#!/usr/bin/env python3
"""Propose nine-wall definitions and rational witnesses for the polar bank."""
from hashlib import sha256
import json
import numpy as np
from scipy.linalg import qr
from fresh_h9_makarov_polar_closure_20261013 import setup, ROOT, TAG, SOURCE, SCOUT


def run():
    source, Q, R, c, refs, points = setup()
    trace_path = ROOT/'results'/f'{TAG}_dfs.json'
    trace = json.loads(trace_path.read_text())
    polar_recipes = []
    for row in refs[27:]:
        ids = row['support']
        matrix = R[ids]@Q
        _, triangular, pivots = qr(matrix.T, mode='economic', pivoting=True)
        assert sum(abs(np.diag(triangular)) > 1e-9) == 9
        support = sorted(ids[int(i)] for i in pivots[:9])
        polar_recipes.append({'scout_vertex': row['polar_vertex'], 'nine_walls': support})
    critical = sorted({row['witness_index'] for row in trace['branch_trace']})
    witnesses = []
    for i in critical:
        point = (1-1e-6)*points[i]+1e-6*c
        integer = list(map(int, np.rint(point*1e12)))
        witnesses.append({'scout_vertex': refs[27+i]['polar_vertex'], 'integer_point': integer})
    out = {'status': 'LITERAL_POLAR_RECIPES_AND_RATIONAL_WITNESS_PROPOSALS',
           'sources': {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest()
                       for p in [SOURCE, SCOUT, trace_path]},
           'polar_recipes': polar_recipes, 'witnesses': witnesses,
           'scope': 'Proposal data only. The exact verifier must prove rank, spacelikeness, all wall signs, and the incompatibility covering obstruction.'}
    (ROOT/'results'/f'{TAG}_recipes.json').write_text(json.dumps(out, indent=2)+'\n')
    print(len(polar_recipes), 'polar recipes;', len(witnesses), 'witnesses')


if __name__ == '__main__':
    run()

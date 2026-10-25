"""SOURCE PROPOSAL ONLY: adaptive prefix identity and synthetic fixtures.

Not a provider or execution record. No repository modules or data are loaded.
Root may copy the pure helper into the final reader after separate review.
Fixtures are definitions only; no main/test runner is supplied or executed.
"""

from copy import deepcopy
from hashlib import sha256
import json
import types
import unittest


def need(ok, message):
    if not ok:
        raise ValueError(message)


def adaptive_prefix_identity(p, core, row, c):
    """Bind an ALREADY F.actual-validated adaptive row to its actual corner.

    Returns identity only. Whole-support geometry and independent same-face
    S8 must remain separate accepted inputs in the enclosing application.
    """
    method = 'nine_context_adaptive_finite_payment'
    key = row['consumer_key']
    fields = {'high_index', 'method', 'core_key', 'local_case', 'selected_marks'}
    need(type(key) is dict and set(key) == fields,
         'complete adaptive consumer key, including local source case')
    need(p['method'] == method and type(p['index']) is int and
         0 <= p['index'] < 2614 and type(p['core_key']) is str and
         type(p['local_case']) is int and 0 <= p['local_case'] < 9,
         'actual adaptive high/core/local-case identity')
    gram = core['gram']
    need(type(gram) is list and 1 <= len(gram) <= 64 and
         all(type(r) is list and len(r) == len(gram) for r in gram),
         'bounded complete physical core Gram')
    marks = key['selected_marks']
    need(type(marks) is list and len(marks) == 2 and
         all(type(v) is int and 0 <= v < len(gram) for v in marks) and
         len(set(marks)) == 2 and type(key['high_index']) is int and
         type(key['local_case']) is int,
         'ordered distinct adaptive actual marks and strict integer key')
    consumer = dict(high_index=p['index'], method=method,
                    core_key=p['core_key'], local_case=p['local_case'],
                    selected_marks=marks)
    need(key == consumer and type(p['payment']) is dict and
         p['payment']['corner'] == marks,
         'same adaptive high/method/core/local-case and actual payment marks')
    corners = p['all_corners']
    need(type(corners) is list and len(corners) == 4 and
         all(type(q) is dict for q in corners),
         'four original actual corner records')
    matches = [(i, q) for i, q in enumerate(corners)
               if q.get('preferred_corner') == marks]
    need(len(matches) == 1, 'unique actual adaptive paying-corner selector')
    ordinal, chosen = matches[0]
    nodes = chosen['marked_component']
    need(type(nodes) is list and len(nodes) == len(set(nodes)) and
         all(type(v) is int and 0 <= v < len(gram) for v in nodes) and
         set(marks) <= set(nodes), 'complete selected adaptive component nodes')
    expected = sha256(c.canonical([gram, nodes, marks]).encode()).hexdigest()
    need(row['physical_context_sha256'] == expected,
         'same adaptive actual core/component/ordered-mark digest')
    local = p['local_case']
    branch_method = ('sealed_marked_leaf' if local in (1, 3, 5, 8) else
                     'adaptive_untouched_marked_leaf' if local in (0, 4, 7) else
                     'adaptive_leaf_or_bounded_relevant_tail')
    need(row['method'] == branch_method and
         type(row['current_alternative_branches']) is list and
         row['current_alternative_branches'] and
         row['required_independent_theorem'] == 'finite57.current.S8.whole_face.v1' and
         row['independent_S8_verified'] is False,
         'same accepted adaptive branch; independent S8 remains a prerequisite')
    path = row['retained_path']
    need(type(path) is list and len(path) == 10 and len(set(path)) == 10 and
         all(type(v) is int and v in nodes for v in path) and
         [path[0], path[-1]] == marks,
         'same ordered ten-node retained source path')
    selector = dict(high_row_selector=['isolation', 'cases', p['index']],
                    actual_corner_selector=['all_corners', ordinal],
                    component_nodes=nodes, selected_marks=marks,
                    local_case=local, identity_only=True,
                    whole_face_payment_proved_here=False,
                    intrinsic_corner_unit=False)
    return consumer, expected, selector


IDENTITY = adaptive_prefix_identity
CANONICAL = types.SimpleNamespace(canonical=lambda value: json.dumps(
    value, sort_keys=True, separators=(',', ':')))


def fixture():
    """Actual adaptive SHAPE only; this is no real core or polytope."""
    gram = [['2' if i == j else '-1' if abs(i-j) == 1 else '0'
             for j in range(12)] for i in range(12)]
    marks = [1, 10]
    nodes = list(range(12))
    p = dict(index=1610, method='nine_context_adaptive_finite_payment',
             core_key='SYNTHETIC_NO_REALIZATION', local_case=5,
             payment=dict(corner=marks), all_corners=[
                 dict(preferred_corner=marks, marked_component=nodes),
                 dict(preferred_corner=[1, 11], marked_component=None),
                 dict(preferred_corner=[2, 10], marked_component=nodes),
                 dict(preferred_corner=[2, 11], marked_component=None)])
    row = dict(consumer_key=dict(high_index=1610, method=p['method'],
               core_key=p['core_key'], local_case=5, selected_marks=marks),
               method='sealed_marked_leaf', retained_path=list(range(1, 11)),
               current_alternative_branches=[dict(basis='SYNTHETIC_CONDITIONAL_ONLY')],
               required_independent_theorem='finite57.current.S8.whole_face.v1',
               independent_S8_verified=False,
               physical_context_sha256=sha256(CANONICAL.canonical(
                   [gram, nodes, marks]).encode()).hexdigest())
    return p, dict(gram=gram), row


class AdaptivePrefixIdentityControls(unittest.TestCase):
    def test_no_invented_top_level_corner_and_exact_actual_selector(self):
        p, core, row = fixture()
        self.assertNotIn('preferred_corner', p)
        self.assertNotIn('marked_component', p)
        consumer, digest, chosen = IDENTITY(p, core, row, CANONICAL)
        self.assertEqual(consumer, row['consumer_key'])
        self.assertEqual(digest, row['physical_context_sha256'])
        self.assertEqual(chosen['actual_corner_selector'], ['all_corners', 0])
        self.assertFalse(chosen['whole_face_payment_proved_here'])
        self.assertFalse(chosen['intrinsic_corner_unit'])

    def test_wrong_core_local_case_or_ordered_marks_refused(self):
        for field, value in (('core_key', 'OTHER_CORE'), ('local_case', 4),
                             ('local_case', True), ('selected_marks', [10, 1])):
            p, core, row = fixture()
            row['consumer_key'][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                IDENTITY(p, core, row, CANONICAL)

    def test_missing_duplicate_or_unmarked_actual_corner_refused(self):
        for mode in ('missing', 'duplicate', 'no_component'):
            p, core, row = fixture()
            if mode == 'missing':
                p['all_corners'][0]['preferred_corner'] = [0, 10]
            elif mode == 'duplicate':
                p['all_corners'][1] = deepcopy(p['all_corners'][0])
            else:
                p['all_corners'][0]['marked_component'] = None
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                IDENTITY(p, core, row, CANONICAL)

    def test_changed_physical_gram_or_digest_refused(self):
        for mode in ('gram', 'digest'):
            p, core, row = fixture()
            if mode == 'gram':
                core['gram'][0][0] = '3'
            else:
                row['physical_context_sha256'] = '0'*64
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                IDENTITY(p, core, row, CANONICAL)

    def test_s8_not_promoted_and_retained_path_keeps_actual_endpoints(self):
        for mode in ('verified', 'wrong_theorem', 'reversed_path', 'wrong_branch'):
            p, core, row = fixture()
            if mode == 'verified':
                row['independent_S8_verified'] = True
            elif mode == 'wrong_theorem':
                row['required_independent_theorem'] = 'UNBOUND_OTHER_PROVIDER'
            elif mode == 'reversed_path':
                row['retained_path'].reverse()
            else:
                row['method'] = 'adaptive_untouched_marked_leaf'
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                IDENTITY(p, core, row, CANONICAL)

"""Require a complete, hash-bound preflight before a research SDPA launch."""
import hashlib
import json
from pathlib import Path


AUDIT_FIELDS={'equalities','pending_face_equations','constant_cones','zero_PSD_dimensions',
              'duplicate_cones','unused_variables'}


def file_digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024**2),b''):h.update(chunk)
    return h.hexdigest()


def scientific_contract(problem):
    m=problem['metadata']
    labels=set(m.get('equality_labels',[]))
    required={f'feature_marginal_phi{i}_potential0' for i in range(26)}
    if m.get('intrinsic_feature_coordinates')!=26 or not required<=labels:
        raise ValueError('The complete 26-feature interface is missing')
    if m.get('selected_vertex_or_face_diagram_lists_used') is not False:
        raise ValueError('Declare complete global families, without selected diagram lists')
    if not m.get('all_measures_are_necessary_moments_not_realizability_assertions'):
        raise ValueError('The necessary-data interpretation is missing')
    if not m.get('actual_polytope_lift') or not m.get('spherical_and_affine_rank_parameters_unbounded'):
        raise ValueError('The uniform actual-polytope lift is missing')
    return dict(intrinsic_feature_coordinates=26,all_feature_marginals_present=True,
        complete_family_scope_declared=True,actual_polytope_lift=m['actual_polytope_lift'],
        threshold=m['threshold'],parity=m['parity'],upper_dimension=m.get('upper_dimension'),
        hierarchy_level=m.get('uniform_strengthening_level',1),
        independent_geometric_review_complete=m.get('independent_source_review_complete',False))


def check_registered_input(input_sha,registration=None):
    root=Path(__file__).resolve().parents[2]
    policy=root/'research/global_formulation/constraint_policy.json'
    rule=json.loads(policy.read_text())
    if input_sha==rule['pinned_research_input_sha256']:
        return dict(kind='pinned_complete_global_model',policy=str(policy.relative_to(root)),
                    policy_sha256=file_digest(policy),bound_artifacts={str(policy.relative_to(root)):file_digest(policy)})
    if registration is None:
        raise ValueError('Changed scientific constraints require a checked constraint-change registration')
    path=Path(registration).resolve();record=json.loads(path.read_text())
    if record.get('schema')!='verified-constraint-change-v1' or record.get('status')!='CONSTRAINT_EXTENSION_CHECKED':
        raise ValueError('Incomplete constraint-change registration')
    if record.get('new_input_sha256')!=input_sha or not record.get('all_parent_constraints_retained'):
        raise ValueError('Registration belongs to a different model or lost old constraints')
    for name,sha in record['bound_artifacts'].items():
        p=(root/name).resolve()
        if not p.is_relative_to(root) or file_digest(p)!=sha:
            raise ValueError('Stale constraint-change source or test receipt: '+name)
    return dict(kind='checked_constraint_extension',registration=str(path.relative_to(root)),
                registration_sha256=file_digest(path),
                bound_artifacts={**record['bound_artifacts'],str(path.relative_to(root)):file_digest(path),
                                 str(policy.relative_to(root)):file_digest(policy)})


def require_clean_audit(audit):
    if audit.get('schema')!='conic-fixed-point-audit-v1' or set(audit.get('counts',{}))!=AUDIT_FIELDS:
        raise ValueError('Missing or incomplete fixed-point audit')
    if audit.get('passed') is not True or any(audit['counts'].values()):
        raise ValueError('Presolve has not reached a verified fixed point')


def native_shape(path):
    with Path(path).open() as stream:
        comment=stream.readline();m=int(stream.readline());n=int(stream.readline())
        sizes=[int(v) for v in stream.readline().split()]
    if not comment.startswith('"Global Coxeter SDP') or len(sizes)!=n or not all(sizes):
        raise ValueError('Unexpected native input shape')
    return m,sizes


def verify_native(input_path,receipt_path,*,fixture=False):
    input_path=Path(input_path).resolve()
    m,sizes=native_shape(input_path)
    if fixture:
        if m>16 or max(abs(r) for r in sizes)>64 or input_path.stat().st_size>1024**2:
            raise ValueError('Only small regression fixtures may bypass the research preflight')
        return dict(scope='small_explicit_regression_fixture',native_input_sha256=file_digest(input_path))
    if receipt_path is None:
        raise ValueError('Research SDPA requires --preflight from the verified preparation pipeline')
    receipt_path=Path(receipt_path).resolve()
    receipt=json.loads(receipt_path.read_text())
    if receipt.get('schema')!='solver-ready-v1' or receipt.get('status')!='SOLVER_READY':
        raise ValueError('The supplied checkpoint is not solver-ready')
    for key in ('presolve_fixed_point','all_exact_reduction_proofs_verified','all_scaling_inverse_identities_verified'):
        if receipt.get(key) is not True:raise ValueError('Missing solver gate: '+key)
    require_clean_audit(receipt['fixed_point_audit'])
    if receipt['native_input_sha256']!=file_digest(input_path):
        raise ValueError('The native matrices changed after preflight')
    if receipt.get('native_free_variables')!=m or receipt.get('native_block_sizes')!=sizes:
        raise ValueError('The preflight describes a different SDP layout')
    root=Path(__file__).resolve().parents[2]
    for name,expected in receipt['bound_artifacts'].items():
        path=(root/name).resolve()
        if not path.is_relative_to(root) or file_digest(path)!=expected:
            raise ValueError('A preflight-bound model, proof or reconstruction changed: '+name)
    if receipt.get('scientific_contract',{}).get('intrinsic_feature_coordinates')!=26:
        raise ValueError('The original feature interface was not checked')
    return dict(scope='verified_global_research_model',preflight=str(receipt_path),
        preflight_sha256=file_digest(receipt_path),native_input_sha256=receipt['native_input_sha256'],
        frozen_input_sha256=receipt['frozen_input_sha256'])

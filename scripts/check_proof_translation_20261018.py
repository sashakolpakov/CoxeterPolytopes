#!/usr/bin/env python3
"""Build and audit the proof-translation package with explicit partial scope.

This does not replay every historical geometry computation and does not
claim end-to-end Lean verification of the geometric research theorems.
"""
import argparse,json,re,sys
from pathlib import Path
from run_lean_batches_signal_safe_20261018 import guarded_run, STATE

ROOT=Path(__file__).resolve().parents[1]
MANUSCRIPTS=('finite_volume58','compact_arithmetic21','ideal_ra4_second_volume',
    'ideal_ra4_second_volume_structure','right_angled_five_face_carriers',
    'right_angled_five_dimensional_constraints','finite_volume57',
    'ideal_right_angled_six_face_constraints','ideal_ra5_cusp_compatibility',
    'compact29_triangle_constraints','ideal_ra5_octahedral_blocks',
    'compact21_supported_triangle_constraints','right_angled_mixed_corner_constraints')

def run(args,cwd=ROOT):
    proof_runner=len(args)>1 and Path(args[1]).name in (
        'run_lean_batches_20261018.py','run_lean_batches_signal_safe_20261018.py')
    try:
        # The runner times each module separately. Give it time to reap its
        # own process group if this outer monitor has to terminate it.
        result=guarded_run(args,cwd=cwd,timeout_seconds=0 if proof_runner else None,
                           termination_grace=10 if proof_runner else 5)
    except (RuntimeError,OSError):
        failure=STATE/'last_failure.txt'
        if failure.exists():
            with failure.open('rb') as stream:
                stream.seek(max(0,failure.stat().st_size-8192))
                print(stream.read(8192).decode(errors='replace'),end='',file=sys.stderr)
        raise
    if result['output_bytes']>len(result['output']):
        raise RuntimeError('Command output exceeded 512KiB; incomplete audit output rejected')
    return result['output'].decode(errors='replace')

def clean_tex_intermediates():
    for pattern in ('*.aux','*.fls','*.fdb_latexmk','*.log','*.out','*.toc','*.synctex.gz'):
        for path in (ROOT/'manuscripts/build').glob(pattern):
            path.unlink()

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--lean-only',action='store_true');args=parser.parse_args()
    registry=json.loads((ROOT/'formal/assumptions.json').read_text())
    assert registry['custom_lean_axioms_declared']==[]
    sources=sorted(p for p in(ROOT/'formal').rglob('*.lean')if'.lake'not in p.parts)
    for path in sources:
        with path.open() as stream:
            for line in stream:
                assert not re.search(r'\b(sorry|admit|native_decide|unsafe)\b',line),path
                assert not re.search(r'^\s*axiom\b',line),path
    for exporter in ('residual58','finite59','ideal59','type_domination59'):
        print(run([sys.executable,f'scripts/export_lean_{exporter}_20261018.py']),end='',flush=True)
    audit=run([sys.executable,'scripts/run_lean_batches_signal_safe_20261018.py','--audit-axioms'])
    print(audit,end='',flush=True)
    assert 'KERNEL_AXIOM_AUDIT' in audit,'Capped axiom audit did not complete'
    allowed=set(registry['allowed_foundational_axioms']);found=set()
    for names in re.findall(r'depends on axioms:\s*\[([^]]*)\]',audit):
        found.update(x.strip()for x in names.split(',')if x.strip())
    assert found,'Expected root theorem axiom reports were missing'
    assert found<=allowed,(found,allowed)
    assert 'sorryAx'not in audit and'Lean.trustCompiler'not in audit
    pdfs=[]
    if not args.lean_only:
        for command in (
            ['scripts/fresh_a15_d15_dimension58_20261018.py','--verify'],
            ['scripts/fresh_d13_context_dimension57_20261018.py','--verify'],
            ['scripts/fresh_arithmetic_dimension21_composition_20261016.py','--bound','21','--verify'],
            ['scripts/fresh_arithmetic_dimension21_scope_check_20261016.py','--verify'],
            ['scripts/check_ideal_ra4_second_volume_saved_20261018.py'],
            ['scripts/check_ideal_ra4_equality_structure_saved_20261018.py'],
            ['scripts/check_ideal_ra4_positive_kappa_saved_20261018.py'],
            ['scripts/fresh_right_angled_mixed_five_face_density_20261018.py','--verify'],
            ['scripts/fresh_right_angled_ideal6_carrier_propagation_20261018.py','--verify'],
            ['scripts/check_ideal_ra5_carrier_cusps_saved_20261018.py'],
            ['scripts/fresh_ideal_ra5_octahedral_blocks_20261018.py','--verify'],
            ['scripts/check_ideal_ra5_octahedral_blocks_saved_20261018.py'],
            ['scripts/check_compact21_e8_triangle_bundles_20261018.py'],
            ['scripts/fresh_compact21_d5_a1_triangle_coupling_20261018.py','--verify'],
            ['scripts/fresh_compact21_supported_forest_triangle_coupling_20261018.py','--verify'],
            ['scripts/check_right_angled_mixed_weighted_triangles_research_20261018.py'],
            ['scripts/check_right_angled_mixed_low_facet_rigidity_research_20261018.py'],
            ['scripts/check_right_angled_mixed_incidence_feasibility_research_20261018.py'],
            ['scripts/check_right_angled_mixed_complete_five_corner_research_20261018.py'],
            ['scripts/check_right_angled_mixed_opposition_moments_research_20261018.py'],
            ['scripts/check_compact29_triangle_star_saved_20261018.py']):
            print(run([sys.executable,*command]),end='',flush=True)
        for name in MANUSCRIPTS:
            source=ROOT/f'manuscripts/{name}.tex'
            assert source.is_file(),f'Manuscript still pending: {source}'
            try:
                run(['latexmk','-pdf','-interaction=nonstopmode','-halt-on-error',
                    '-outdir=manuscripts/build',str(source.relative_to(ROOT))])
            finally:
                clean_tex_intermediates()
            pdf=ROOT/f'manuscripts/build/{name}.pdf';assert pdf.is_file()and pdf.stat().st_size>0
            pdfs.append(str(pdf.relative_to(ROOT)))
    print(json.dumps(dict(status='PASS_PARTIAL_PROOF_TRANSLATION',source_files_audited=len(sources),
        foundational_axioms=sorted(found),custom_axioms=[],unfinished_proof_terms=False,
        real_certificates=['149-item rank58 residual upper378',
            'all5148 ambient59 residual tables at used caps',
            'all1890 terminal ideal profiles and complete terminal list',
            '180 actual/representative bases and32400 cross pairs'],pdfs=pdfs,
        full_geometric_theorems_lean_verified=False),indent=2))

if __name__=='__main__':main()

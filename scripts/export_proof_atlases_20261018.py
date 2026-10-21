#!/usr/bin/env python3
"""Translate the bound proof-report dependency archives to TeX via Pandoc.

This is an archival translation, not a replacement for mathematical review.
The frozen files and their historical status statements are preserved.
"""
import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PACKAGES={
 'finite_volume57':('results/fresh_d13_context_dimension57_20261018.json',
     'docs/finite_volume57_lemma_map.md'),
 'finite_volume58':('results/fresh_a15_d15_dimension58_20261018.json',
     'docs/finite_volume58_lemma_map.md'),
 'compact_arithmetic21':('results/fresh_arithmetic_dimension21_composition_20261016_bound21.json',
     'docs/compact_arithmetic21_proof_map.md'),
 'right_angled':([
     'results/fresh_ideal_ra4_double24_volume_20261018.json',
     'results/fresh_ideal_ra4_equality_structure_20261018.json',
     'results/fresh_right_angled_five_face_carrier_20261018.json',
     'results/fresh_ideal_ra4_positive_kappa_doubling_20261018.json',
     'results/fresh_ideal_ra5_positive_carrier_one_layer_20261018.json',
     'results/fresh_right_angled_mixed_five_face_density_20261018.json',
     'results/fresh_right_angled_ideal6_carrier_propagation_20261018.json',
     'results/fresh_ideal_ra5_known_carrier_cusp_compatibility_20261018.json',
     'results/fresh_ideal_ra5_octahedral_blocks_20261018.json'],
     'docs/right_angled_proof_index.md'),
 'compact_research':([
     'results/fresh_compact29_three_triangle_star_constraint_20261018.json',
     'results/fresh_compact21_triangle_transfer_scope_20261018.json',
     'results/fresh_exact66_old7_boundary_cap_obstruction_20261018.json'],
     'docs/compact_research_proof_index.md')}

def digest(p):
    with p.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def streamed_manifest_fields(path):
    # jq --stream retains paths/scalars, not the potentially huge data arrays.
    # Only top-level source/hash-map entries are needed to walk the archive.
    query='''reduce inputs as $i ({sources:{},source_hashes:{}};
      if (($i|length)==2 and ($i[0]|length)==2 and
          ($i[0][0]=="sources" or $i[0][0]=="source_hashes"))
      then .[$i[0][0]][$i[0][1]]=$i[1] else . end)'''
    completed=subprocess.run(['jq','--stream','-n',query,str(path)],
        capture_output=True,text=True,check=True,timeout=60)
    if len(completed.stdout)>524288:
        raise RuntimeError(f'{path}: source manifest exceeds512KiB')
    return json.loads(completed.stdout)

def collect(endpoint,proofmap):
    seen=set();reports=set();stack=[endpoint]if isinstance(endpoint,str)else list(endpoint)
    while stack:
        name=stack.pop()
        if name in seen:continue
        seen.add(name);path=ROOT/name
        if name.endswith('.md'):reports.add(name)
        if not name.endswith('.json'):continue
        if path.stat().st_size>20*1024*1024:
            record=streamed_manifest_fields(path)
        else:
            with path.open() as stream:
                record=json.load(stream)
        for field in ('sources','source_hashes'):
            for source,expected in record.get(field,{}).items():
                assert digest(ROOT/source)==expected,(name,source)
                stack.append(source)
    # The proof map also names prose inputs not stored in an old JSON manifest.
    for link in re.findall(r'\]\(([^)]+\.md)\)',(ROOT/proofmap).read_text()):
        path=((ROOT/proofmap).parent/link).resolve()
        if path.is_file()and path.is_relative_to(ROOT):reports.add(str(path.relative_to(ROOT)))
    return sorted(reports),sorted(seen)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    outdir=ROOT/'manuscripts/generated';outdir.mkdir(exist_ok=True)
    for package,(endpoint,proofmap) in PACKAGES.items():
        reports,dependencies=collect(endpoint,proofmap)
        endpoint_display=('`'+endpoint+'`'if isinstance(endpoint,str)else
            '\n'.join('* `'+p+'`'for p in endpoint))
        intro=f'''# Frozen proof-report archive: {package}

This is a mechanical TeX translation of every bound Markdown proof report,
plus prose reports explicitly named by the current proof map. It accompanies
the edited manuscript, exact scripts and partial Lean development. It is
an archival source collection, not a claim of end-to-end formal verification.

Reports retain their original dates, scopes and historical status wording.
Later reports can supersede earlier bounds. The authoritative endpoint inputs are:

{endpoint_display}

The current edited proof map is `{proofmap}`.

The manifest records all report hashes and every recursively bound file.
This export verifies source identity; it does not reprove the geometry.

'''
        pieces=[intro]
        for i,name in enumerate(reports,1):
            body=(ROOT/name).read_text()
            # Prefix the source path; keep the complete original report below it.
            pieces.append(f'\n\\clearpage\n\n# Source {i}\n\n'
                f'Source path: \\path{{{name}}}\n\n'
                f'SHA256: `{digest(ROOT/name)}`\n\n'+body+'\n')
        parsed=json.loads(subprocess.run(['pandoc',
            '--from=markdown-yaml_metadata_block+tex_math_single_backslash','--to=json'],
            input='\n'.join(pieces),text=True,capture_output=True,check=True).stdout)
        def protect_bare_math(node):
            if isinstance(node,dict):
                if node.get('t')=='Math'and re.search(r'\\[()]',node['c'][1]):
                    raw=node['c'][1]
                    node.update(t='Code',c=[['',[],[]],raw])
                if node.get('t')=='RawInline'and node['c'][0]=='tex':
                    raw=node['c'][1]
                    command=re.match(r'\\([A-Za-z]+)',raw)
                    text_commands={
                        'path','url','nolinkurl','href','ref','label','texttt',
                        'textbf','textit','emph','text','quad','qquad','newline'}
                    math_commands=set(('binom frac sqrt sum prod ker det rank mathbb mathcal '
                        'langle rangle in notin le ge leq geq ne neq pm mp cdot times to mapsto '
                        'longrightarrow Longrightarrow iff implies subset subseteq supset supseteq '
                        'cup cap setminus smallsetminus sigma alpha beta gamma delta epsilon '
                        'varepsilon zeta eta theta vartheta iota kappa lambda mu nu xi omicron '
                        'pi rho varrho tau upsilon phi varphi chi psi omega Gamma Delta Theta '
                        'Lambda Xi Pi Sigma Upsilon Phi Psi Omega tilde widehat hat overline '
                        'underline left right begin end operatorname mathrm mathbf mathsf '
                        'mathtt textstyle displaystyle limits nolimits tag ensuremath '
                        'ell infinity infty partial nabla vert Vert lvert rvert lVert rVert '
                        'lfloor rfloor lceil rceil').split())
                    if command and command.group(1)not in text_commands:
                        needs_argument={'binom','frac','sqrt','mathbb','mathcal','mathsf',
                            'mathbf','mathtt','mathrm','operatorname','tilde','widehat',
                            'hat','overline','underline','left','right','tag'}
                        bare_argument_macro=(command.group(1)in needs_argument and
                            raw.strip()=='\\'+command.group(1))
                        if command.group(1)in math_commands and not bare_argument_macro:
                            node['c'][1]='\\ensuremath{'+raw+'}'
                        else:
                            # An informal backslash such as F\\a is source
                            # text, not an invented LaTeX macro definition.
                            node.update(t='Code',c=[['',[],[]],raw])
                if node.get('t')=='Code'and len(node['c'][1])>40 and re.fullmatch(
                        r'[A-Za-z0-9_./+\-]+',node['c'][1]):
                    node.update(t='RawInline',c=['tex','\\nolinkurl{'+node['c'][1]+'}'])
                for value in node.values():protect_bare_math(value)
            elif isinstance(node,list):
                for value in node:protect_bare_math(value)
        protect_bare_math(parsed)
        converted=subprocess.run(['pandoc','--from=json',
            '--to=latex','--standalone','--top-level-division=chapter',
            '--variable=documentclass:report','--variable=mainfont:STIXGeneral',
            '--variable=mathfont:STIX Two Math','--variable=monofont:Menlo',
            '--variable=geometry:margin=22mm','--variable=fontsize:10pt',
            '--metadata=title:'+package+' frozen proof archive'],input=json.dumps(parsed),
            text=True,capture_output=True,check=True).stdout
        # Preserve the exact Unicode symbols absent from STIXGeneral using
        # a font which contains them. Source mathematics is not rewritten.
        fallback='\\newfontfamily\\ArchiveFallback{Arial Unicode MS}\n'
        fallback+='\\usepackage{newunicodechar}\n'
        for char in '⁻₀₁₂₃₄₅₆₇₈₉₋₊':
            fallback+=f'\\newunicodechar{{{char}}}{{{{\\ArchiveFallback {char}}}}}\n'
        for char,letter in zip('ₙᵢᵣᵥᵤₐₑₒₓₕₖₗₘₚₛₜ','nirvuaeoxhklmpst'):
            fallback+='\\newunicodechar{'+char+'}{\\ensuremath{{}_{\\mathrm{'+letter+'}}}}\n'
        fallback+='\\newunicodechar{ᵀ}{\\ensuremath{{}^{\\mathsf{T}}}}\n'
        fallback+='\\newunicodechar{ℓ}{\\ensuremath{\\ell}}\n'
        fallback+='\\usepackage{xurl}\n\\usepackage{fvextra}\n\\fvset{breaklines=true,breakanywhere=true}\n'
        fallback+='\\hypersetup{hypertexnames=false}\n'
        converted=converted.replace('\\begin{document}',fallback+'\\begin{document}',1)
        tex=outdir/f'{package}_atlas.tex'
        manifest=outdir/f'{package}_atlas_manifest.json'
        data={'status':'MECHANICAL_TEX_ARCHIVE_NOT_NEW_PROOF',
            'endpoint':endpoint,'proof_map':proofmap,
            'reports':{name:digest(ROOT/name) for name in reports},
            'recursive_dependencies':{name:digest(ROOT/name) for name in dependencies},
            'tex_sha256':hashlib.sha256(converted.encode()).hexdigest(),
            'pandoc_version':subprocess.run(['pandoc','--version'],text=True,capture_output=True,
                check=True).stdout.splitlines()[0]}
        if args.verify:
            assert tex.read_text()==converted
            assert json.loads(manifest.read_text())==data
        else:
            tex.write_text(converted);manifest.write_text(json.dumps(data,indent=2)+'\n')
        print(package,len(reports),'complete reports,',len(dependencies),'bound files',flush=True)

if __name__=='__main__':main()

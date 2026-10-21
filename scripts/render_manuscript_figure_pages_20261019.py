#!/usr/bin/env python3
"""Render bounded manuscript figure pages for visual review, not proof.

Run only inside the coordinator's shared-lock/resource-capped wrapper.
The output PNGs are disposable; retain a visual-review ledger and remove
them after inspection. The input PDF is never modified.
"""
import argparse
from hashlib import file_digest
import json
from pathlib import Path
import re
import signal
import subprocess
import tempfile
from pypdf import PdfReader, PdfWriter
from pypdf.generic import DecodedStreamObject

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pdf')
    parser.add_argument('--sha256', required=True)
    parser.add_argument('--output-name', required=True)
    args = parser.parse_args()
    signal.alarm(45)
    pdf = (ROOT/args.pdf).resolve()
    assert pdf.is_relative_to(ROOT/'manuscripts/build')
    assert 0 < pdf.stat().st_size < 4*1024**2
    assert re.fullmatch(r'[a-z0-9_]{1,80}', args.output_name)
    out = ROOT/'formal/.proof-batches'/args.output_name
    assert not out.exists()
    with pdf.open('rb') as stream:
        assert file_digest(stream, 'sha256').hexdigest() == args.sha256
    reader = PdfReader(pdf)
    assert 1 <= len(reader.pages) <= 24
    pages = []
    for i, page in enumerate(reader.pages):
        figures = re.findall(r'Figure\s+(\d+)\s*:', page.extract_text())
        if figures:
            pages.append((i, list(map(int, figures))))
    assert 1 <= len(pages) <= 8
    out.mkdir()
    rows = []
    with tempfile.TemporaryDirectory(prefix='manuscript-figures-', dir='/private/tmp') as temp:
        for i, figures in pages:
            writer = PdfWriter()
            writer.add_page(reader.pages[i])
            canvas = writer.pages[0]
            stream = DecodedStreamObject()
            width, height = float(canvas.mediabox.width), float(canvas.mediabox.height)
            stream.set_data(f'q 1 1 1 rg 0 0 {width} {height} re f Q\n'.encode()
                            + canvas.get_contents().get_data())
            canvas.replace_contents(stream)
            single = Path(temp)/f'page{i+1}.pdf'
            writer.write(single)
            target = out/f'page{i+1}.png'
            result = subprocess.run(['/usr/bin/sips', '-s', 'format', 'png',
                                     '-Z', '1400', str(single), '--out', str(target)],
                                    capture_output=True, timeout=8)
            assert result.returncode == 0, result.stderr[-1000:]
            assert 0 < target.stat().st_size < 2*1024**2
            rows.append(dict(figures=figures, page=i+1, bytes=target.stat().st_size,
                             png=str(target.relative_to(ROOT))))
    signal.alarm(0)
    print(json.dumps(dict(status='RENDERED_FOR_VISUAL_REVIEW_NOT_PROOF',
                         pdf_sha256=args.sha256, figures=rows)))


if __name__ == '__main__':
    main()

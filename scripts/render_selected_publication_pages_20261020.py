#!/usr/bin/env python3
"""Render selected PDF pages for visual inspection, never mathematical proof.

Use the existing bounded research wrapper with --parent-lock. Delete the
output PNGs after review; the input PDF is read only.
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
    parser.add_argument('--pages', required=True, help='One-based comma-separated pages')
    parser.add_argument('--output-name', required=True)
    args = parser.parse_args()
    signal.alarm(45)
    pdf = (ROOT / args.pdf).resolve()
    assert pdf.is_relative_to(ROOT / 'manuscripts/build')
    assert 0 < pdf.stat().st_size < 4 * 1024**2
    assert re.fullmatch(r'[a-z0-9_]{1,80}', args.output_name)
    assert re.fullmatch(r'\d+(,\d+){0,5}', args.pages)
    pages = [int(p) for p in args.pages.split(',')]
    assert len(set(pages)) == len(pages)
    with pdf.open('rb') as stream:
        assert file_digest(stream, 'sha256').hexdigest() == args.sha256
    reader = PdfReader(pdf)
    # The compact working articles include long case appendices. Keep the
    # six-page render batch unchanged while allowing selection within them.
    assert 1 <= len(reader.pages) <= 512
    assert all(1 <= p <= len(reader.pages) for p in pages)
    out = ROOT / 'formal/.proof-batches' / args.output_name
    out.mkdir(exist_ok=False)
    rows = []
    with tempfile.TemporaryDirectory(prefix='publication-pages-', dir='/private/tmp') as temp:
        for number in pages:
            writer = PdfWriter()
            writer.add_page(reader.pages[number - 1])
            canvas = writer.pages[0]
            stream = DecodedStreamObject()
            width, height = float(canvas.mediabox.width), float(canvas.mediabox.height)
            stream.set_data(f'q 1 1 1 rg 0 0 {width} {height} re f Q\n'.encode()
                            + canvas.get_contents().get_data())
            canvas.replace_contents(stream)
            single = Path(temp) / f'page{number}.pdf'
            writer.write(single)
            target = out / f'page{number}.png'
            result = subprocess.run(['/usr/bin/sips', '-s', 'format', 'png',
                                     '-Z', '1400', str(single), '--out', str(target)],
                                    capture_output=True, timeout=8)
            assert result.returncode == 0, result.stderr[-1000:]
            assert 0 < target.stat().st_size < 2 * 1024**2
            rows.append({'page': number, 'png': str(target.relative_to(ROOT)),
                         'bytes': target.stat().st_size})
    signal.alarm(0)
    print(json.dumps({'status': 'RENDERED_FOR_VISUAL_REVIEW_NOT_PROOF',
                      'pdf_sha256': args.sha256, 'pages': rows}))


if __name__ == '__main__':
    main()

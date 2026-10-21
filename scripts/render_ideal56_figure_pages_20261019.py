#!/usr/bin/env python3
"""Render the four figure pages for manual review; coordinator wrapper only.

This is typesetting inspection, not mathematical verification. The four
small PNG files are disposable after inspection; the PDF remains unchanged.
"""
from hashlib import file_digest
from pathlib import Path
import json
import signal
import subprocess
import tempfile
from pypdf import PdfReader, PdfWriter
from pypdf.generic import DecodedStreamObject

ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT/'manuscripts/build/ideal_right_angled_six_face_constraints.pdf'
EXPECTED = 'a89a13677f3d7256afe5c2837552a669c5dc11e5967119324c6881b7a582d376'
OUT = ROOT/'formal/.proof-batches/ideal56_figure_review_20261019'


def main():
    signal.alarm(40)
    with PDF.open('rb') as stream:
        assert file_digest(stream, 'sha256').hexdigest() == EXPECTED
    reader = PdfReader(PDF)
    assert len(reader.pages) == 12
    assert not OUT.exists()
    OUT.mkdir()
    rows = []
    with tempfile.TemporaryDirectory(prefix='ideal56-figures-', dir='/private/tmp') as temp:
        for figure, page in enumerate((3, 4, 7, 10), 1):
            assert f'Figure {figure}:' in reader.pages[page-1].extract_text()
            writer = PdfWriter()
            writer.add_page(reader.pages[page-1])
            # Paint the temporary PDF canvas white before its original vector
            # content. This avoids transparent PNGs in dark image viewers.
            canvas = writer.pages[0]
            stream = DecodedStreamObject()
            width, height = float(canvas.mediabox.width), float(canvas.mediabox.height)
            background = f'q 1 1 1 rg 0 0 {width} {height} re f Q\n'.encode()
            stream.set_data(background + canvas.get_contents().get_data())
            canvas.replace_contents(stream)
            single = Path(temp)/f'page{page}.pdf'
            writer.write(single)
            target = OUT/f'figure{figure}_page{page}.png'
            result = subprocess.run(['/usr/bin/sips', '-s', 'format', 'png',
                '-Z', '1400', str(single), '--out', str(target)],
                capture_output=True, timeout=8)
            assert result.returncode == 0, result.stderr[-1000:]
            assert target.is_file() and target.stat().st_size < 2*1024**2
            rows.append(dict(figure=figure, page=page, png=str(target.relative_to(ROOT)),
                             bytes=target.stat().st_size))
    signal.alarm(0)
    print(json.dumps(dict(status='RENDERED_FOR_MANUAL_REVIEW_NOT_PROOF',
                         pdf_sha256=EXPECTED, figures=rows)))


if __name__ == '__main__':
    main()

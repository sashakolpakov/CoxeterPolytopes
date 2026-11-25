#!/usr/bin/env python3
"""Compress named research inputs, checking the complete byte round trip.

The original inputs stay in place. This is not a bulk upload or a
mathematical verification; it creates small, restorable Git artifacts.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil

ROOT=Path(__file__).resolve().parents[1]


def digest(stream):
    h=hashlib.sha256();size=0
    while chunk:=stream.read(1024**2): h.update(chunk);size+=len(chunk)
    return size,h.hexdigest()


def main():
    if not os.environ.get("COXETER_CHECK_OWNER"): raise RuntimeError("Use the serial scientific watchdog")
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True)
    parser.add_argument("inputs",nargs="+")
    args=parser.parse_args();out=(ROOT/args.output).resolve()
    if not out.is_relative_to(ROOT/"results") or out.exists(): raise ValueError("Use a fresh receipt")
    records=[]
    for name in args.inputs:
        source=(ROOT/name).resolve();target=source.with_suffix(source.suffix+".gz")
        if not source.is_relative_to(ROOT/"results") or source.suffix!=".json" or target.exists(): raise ValueError("Use fresh JSON archive names in results/")
        with source.open("rb") as stream: before=digest(stream)
        with source.open("rb") as stream,target.open("xb") as raw:
            with gzip.GzipFile(filename="",mode="wb",fileobj=raw,compresslevel=9,mtime=0) as compressed:
                shutil.copyfileobj(stream,compressed,length=1024**2)
        with gzip.open(target,"rb") as stream: restored=digest(stream)
        with source.open("rb") as stream: after=digest(stream)
        if before!=restored or before!=after: raise ValueError("Research input byte round trip failed")
        with target.open("rb") as stream: packed=digest(stream)
        records.append(dict(source=str(source.relative_to(ROOT)),archive=str(target.relative_to(ROOT)),
            source_bytes=before[0],source_sha256=before[1],archive_bytes=packed[0],archive_sha256=packed[1],
            complete_byte_round_trip_verified=True,original_preserved=True))
    result=dict(status="RESEARCH_MODEL_BYTE_ROUND_TRIPS_PASSED",records=records,mathematical_verification_asserted=False,upload_performed=False)
    with out.open("x") as stream: json.dump(result,stream,indent=2,sort_keys=True);stream.write("\n")
    print(json.dumps(result),flush=True)


if __name__=="__main__": main()

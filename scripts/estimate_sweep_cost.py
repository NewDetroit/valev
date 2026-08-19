#!/usr/bin/env python3
"""Estimate triage-sweep transfer volume from the Logan per-accession stats parquet.

Reads  : $WORK/accessions.txt          (one SRA run accession per line)
Writes : stdout, and optionally a two-column TSV via --out

Run this BEFORE scripts/05_logan_triage.sh. Two numbers matter:

  * matched fraction — how much of the niche exists in this Logan release at all.
    Accessions absent from Logan are invisible to the fast path, full stop.
  * total bytes      — transfer volume, NOT disk requirement. The sweep streams
    contigs through DIAMOND and never stores them.

CORRECTION C-07 — PLAN.md Task 6 Step 5 selected the size column by substring match:

    col = [c for c in sub.columns if "contig" in c and "after" in c]
    if col: tb = sub[col[0]].sum()/1e12

Two defects. On the v1 parquet that selector matches TWO columns
(`contigs_after_compression` and `size_contigs_after_compression`) and `col[0]` takes
whichever pandas happens to order first, so the figure that decides whether to spend the
compute budget is arbitrary. And `if col:` means a schema change produces no output at all
rather than an error — a silent skip on the one step whose entire job is to stop you
overspending. Here the exact column comes from $LOGAN_SIZE_COL and its absence is fatal.

Verified against both parquet footers on 2026-08-19:
    logan-seqstats-contigs-v1.2.parquet  38,124,741 rows  contigs_after_compression_bytes
    logan-seqstats.parquet (v1)          27,269,310 rows  size_contigs_after_compression

Nulls are reported rather than swallowed: v1.2 covers 38.1M accessions but only 37,377,661
have contigs, and a bare .sum() over the nulls understates the total without saying so.
"""

from __future__ import annotations

import argparse
import io
import os
import sys
import urllib.request

import pyarrow.parquet as pq


def env(name: str) -> str:
    v = os.environ.get(name)
    if not v:
        sys.exit(f"[estimate_sweep_cost] {name} unset — run: source config/config.sh")
    return v


class HTTPRangeFile(io.RawIOBase):
    """Seekable read-only file over HTTP Range requests.

    Parquet keeps its schema and row-group index in a footer and stores columns in
    separate chunks, so a reader that can seek fetches only the footer plus the two
    columns we project. That is what lets a 1.6 GB file be queried without downloading it.
    """

    def __init__(self, url: str, timeout: int = 120):
        self.url, self.timeout, self._pos = url, timeout, 0
        self.n_requests = self.n_bytes = 0
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=timeout) as r:
            if r.headers.get("Accept-Ranges") != "bytes":
                sys.exit(f"[estimate_sweep_cost] {url} does not advertise byte ranges")
            self.size = int(r.headers["Content-Length"])

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self._pos

    def seek(self, offset: int, whence: int = os.SEEK_SET) -> int:
        if whence == os.SEEK_SET:
            self._pos = offset
        elif whence == os.SEEK_CUR:
            self._pos += offset
        elif whence == os.SEEK_END:
            self._pos = self.size + offset
        else:
            raise ValueError(f"bad whence: {whence}")
        return self._pos

    def read(self, size: int = -1) -> bytes:
        if size is None or size < 0:
            size = self.size - self._pos
        size = min(size, self.size - self._pos)
        if size <= 0:
            return b""
        first, last = self._pos, self._pos + size - 1
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={first}-{last}"})
        with urllib.request.urlopen(req, timeout=self.timeout) as r:
            if r.status != 206:
                raise RuntimeError(f"{self.url}: expected HTTP 206, got {r.status}")
            data = r.read()
        self.n_requests += 1
        self.n_bytes += len(data)
        self._pos += len(data)
        return data

    def readinto(self, buf) -> int:
        data = self.read(len(buf))
        buf[: len(data)] = data
        return len(data)


def load_accessions(path: str) -> set[str]:
    if not os.path.exists(path):
        sys.exit(
            f"[estimate_sweep_cost] missing accession list: {path}\n"
            f"  run scripts/04_select_accessions.sh (needs NCBI) or\n"
            f"      scripts/04b_select_accessions_offline.sh (no NCBI)"
        )
    accs = {ln.strip() for ln in open(path) if ln.strip()}
    if not accs:
        sys.exit(f"[estimate_sweep_cost] {path} is empty")
    return accs


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB", "PB"):
        if abs(n) < 1000:
            return f"{n:,.1f} {unit}"
        n /= 1000
    return f"{n:,.1f} EB"


def main() -> int:
    release = env("LOGAN_RELEASE")
    size_col = env("LOGAN_SIZE_COL")
    parquet_name = env("LOGAN_STATS_PARQUET")

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--accessions", default=os.path.join(env("WORK"), "accessions.txt"))
    ap.add_argument("--parquet", help="local stats parquet (skips the network)")
    ap.add_argument(
        "--parquet-url",
        default=f"{env('LOGAN_HTTPS')}/stats/{parquet_name}",
        help="read over HTTP Range instead of downloading (default: the release's parquet)",
    )
    ap.add_argument("--out", help="also write the summary as a two-column TSV")
    args = ap.parse_args()

    accs = load_accessions(args.accessions)

    local = args.parquet or os.path.join(env("REF"), parquet_name)
    handle = None
    if os.path.exists(local):
        source, pf = local, pq.ParquetFile(local)
    else:
        source = args.parquet_url
        handle = HTTPRangeFile(source)
        pf = pq.ParquetFile(handle)

    names = pf.schema_arrow.names
    if size_col not in names:
        sys.stderr.write(
            f"[estimate_sweep_cost] FATAL: column '{size_col}' absent from {source}\n"
            f"  columns present: {', '.join(names)}\n"
            f"  LOGAN_RELEASE={release} expects that column. Either the release moved or\n"
            f"  LOGAN_SIZE_COL in config/config.sh is wrong. This is C-07 — it must fail\n"
            f"  loudly rather than silently reporting nothing.\n"
        )
        return 2

    matched = 0
    nulls = 0
    total = 0
    scanned = 0
    for batch in pf.iter_batches(batch_size=131072, columns=["accession", size_col]):
        acc_col = batch.column(0).to_pylist()
        sz_col = batch.column(1).to_pylist()
        scanned += len(acc_col)
        for a, s in zip(acc_col, sz_col):
            if a in accs:
                matched += 1
                if s is None:
                    nulls += 1
                else:
                    total += int(s)
        if sys.stderr.isatty():
            sys.stderr.write(f"\r  scanned {scanned:,} rows, matched {matched:,}")
            sys.stderr.flush()
    if sys.stderr.isatty():
        sys.stderr.write("\n")

    requested = len(accs)
    frac = matched / requested if requested else 0.0
    sized = matched - nulls

    print()
    print(f"Logan release        : {release}  ({parquet_name})")
    print(f"source               : {source}")
    if handle:
        print(f"  fetched over HTTP  : {handle.n_requests} range requests, "
              f"{human(handle.n_bytes)} (parquet is {human(handle.size)})")
    print(f"parquet rows scanned : {scanned:,}")
    print()
    print(f"accessions requested : {requested:,}")
    print(f"  present in Logan   : {matched:,}  ({frac:.1%})")
    print(f"  absent from Logan  : {requested - matched:,}  "
          f"— invisible to the fast path in this release")
    print(f"  present, no contigs: {nulls:,}  (null {size_col})")
    print(f"  sized              : {sized:,}")
    print()
    print(f"compressed contigs   : {total:,} bytes  = {total/1e12:.3f} TB "
          f"({total/1e9:.1f} GB)")
    print("                       transfer volume, not disk — the sweep streams.")

    if matched and sized:
        mean = total / sized
        print(f"mean per accession   : {human(mean)}")
        if nulls:
            print(f"  NOTE: {nulls:,} matched accessions have no contigs; the total above")
            print( "        covers only the sized ones. Do not present it as a full-set figure.")
    if frac < 0.5 and requested:
        print()
        print(f"WARNING: only {frac:.1%} of the niche is in Logan {release}.")
        print("  Most of this niche postdates the freeze or was never assembled.")
        print("  Consider whether the fast path can answer your question at all.")

    if args.out:
        with open(args.out, "w") as fh:
            fh.write("metric\tvalue\n")
            for k, v in [
                ("logan_release", release), ("parquet", parquet_name),
                ("accessions_requested", requested), ("accessions_matched", matched),
                ("matched_fraction", f"{frac:.6f}"), ("accessions_null_contigs", nulls),
                ("accessions_sized", sized), ("compressed_bytes", total),
                ("compressed_tb", f"{total/1e12:.6f}"),
            ]:
                fh.write(f"{k}\t{v}\n")
        print(f"\nwrote {args.out}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

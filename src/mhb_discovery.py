#!/usr/bin/env python
"""Discovery of methylation haplotype blocks (MHBs) from mHap files.

Re-implementation of MHBDiscovery from mHapSuite (https://github.com/yoyoong/mHapSuite,
com/MHBDiscovery.java and com/common/Util.java), replacing the original mHapTk code.

Regions are scanned in chunks of 1 Mb that overlap by 1 kb, and only CpG sites inside a
chunk are considered. Within a chunk, CpG sites are visited from left to right: site e
extends the current block if it is linked to each of the window - 1 sites before it,
i.e. for every pair (e - i, e), i = 1 .. window - 1, r2 >= r2 cutoff and p <= p cutoff.
When pair (e - i, e) fails, the block [start, e - 1] is reported if it holds at least
`window` sites, and the next block starts after the failing site, at max(start, e - i + 1).
The block left open at the end of the chunk is reported if end - start >= window.
Identical blocks found in several chunks or regions are reported once.

For two CpG sites A < B, r2 is the signed squared correlation of their methylation calls
over the mHaps covering both, and p = P(X >= N11) with X ~ Binomial(N, pA * pB).

Output: one MHB per line, "chrom<TAB>start<TAB>end", where start and end are the 1-based
positions of the first and last CpG site of the block (as in mHapSuite).
"""
import argparse
import gzip
import os
import sys
import time

import numpy as np
from scipy.stats import binom

# mHapSuite processes regions in 1 Mb chunks that overlap by 1 kb
CHUNK_SIZE = 1000000
CHUNK_OVERLAP = 1000


class TabixFile:
    """Tabix reader backed by pysam, or by pytabix (the binding used by the rest of mHapTk)."""

    def __init__(self, path):
        try:
            import pysam
        except ImportError:
            import tabix
            self._tb = tabix.open(path)
            self._error = tabix.TabixError
            self._pysam = False
        else:
            self._tb = pysam.TabixFile(path)
            self._contigs = set(self._tb.contigs)
            self._pysam = True

    def fetch(self, chrom, beg, end):
        """Records overlapping the 0-based half-open interval [beg, end), as lists of fields."""
        if self._pysam:
            if chrom not in self._contigs:
                return
            for line in self._tb.fetch(chrom, beg, end):
                yield line.split('\t')
        else:
            try:
                yield from self._tb.query(chrom, beg, end)
            except self._error:  # chromosome absent from the index
                return


def parse_region(text):
    """'chr:start-end' -> (chrom, start, end), 1-based and closed."""
    chrom, span = text.rsplit(':', 1)
    start, end = span.replace(',', '').split('-')
    return chrom, int(start), int(end)


def read_bed(path):
    """Regions of a BED file, converted to 1-based closed coordinates."""
    regions = []
    with open(path) as fh:
        for line in fh:
            fields = line.split()
            if not fields or fields[0].startswith(('#', 'track', 'browser')):
                continue
            if len(fields) < 3:
                raise ValueError('Interval not in correct format: ' + line.rstrip())
            regions.append((fields[0], int(fields[1]) + 1, int(fields[2])))
    return regions


def genome_regions(cpg_path):
    """One region per chromosome of the CpG file, from its first to its last CpG site."""
    regions = []
    chrom = first = last = None
    with gzip.open(cpg_path, 'rt') as fh:
        for line in fh:
            c, pos = line.split('\t', 2)[:2]
            if c != chrom:
                if chrom is not None:
                    regions.append((chrom, first, int(last)))
                chrom, first = c, int(pos)
            last = pos
    if chrom is not None:
        regions.append((chrom, first, int(last)))
    return regions


def split_region(start, end):
    """Chunks (start, end) of a region, as in mHapSuite's Util.splitRegionToSmallRegion."""
    if end - start <= CHUNK_SIZE:
        return [(start, end)]
    chunks = []
    for _ in range((end - start) // CHUNK_SIZE + 1):
        chunk_end = min(start + CHUNK_SIZE + CHUNK_OVERLAP - 1, end)
        chunks.append((start, chunk_end))
        start = chunk_end - CHUNK_OVERLAP + 1 if chunk_end - CHUNK_OVERLAP + 1 >= 1 else chunk_end + 1
    return chunks


def fetch_mhaps(mhap, chrom, start, end):
    """mHaps overlapping [start, end]; mHaps with missing calls ('.') are skipped as in mHapSuite."""
    starts, ends, haps, counts = [], [], [], []
    for rec in mhap.fetch(chrom, start - 1, end):
        hap = rec[3]
        if '.' in hap:
            continue
        starts.append(int(rec[1]))
        ends.append(int(rec[2]))
        haps.append(hap)
        counts.append(int(rec[4]))
    return np.array(starts, dtype=np.int64), np.array(ends, dtype=np.int64), haps, np.array(counts, dtype=np.int64)


def fetch_cpgs(cpg, chrom, start, end):
    """Sorted positions of the CpG sites in [start, end]."""
    pos = np.array([int(rec[1]) for rec in cpg.fetch(chrom, max(start - 1, 0), end)], dtype=np.int64)
    return pos[(pos >= start) & (pos <= end)]


def pair_counts(first, haps, lens, counts, n_sites, max_dist):
    """Pattern counts of CpG pairs (k, k + d) for d = 1 .. max_dist.

    first[r] is the index of the first CpG site of mHap r. Returns {d: array (n_sites, 4)}
    whose row k holds the counts of the patterns 00, 01, 10 and 11 at (k, k + d).
    """
    meth = np.frombuffer(''.join(haps).encode('ascii'), dtype=np.uint8) - ord('0')
    if meth.size and meth.max() > 1:
        raise ValueError('mHap haplotype strings must consist of 0, 1 and .')
    read = np.repeat(np.arange(len(haps)), lens)
    offset = np.arange(meth.size) - np.repeat(np.cumsum(lens) - lens, lens)
    site = first[read] + offset
    rest = lens[read] - 1 - offset  # sites that follow in the same mHap
    weight = counts[read]
    tables = {}
    for d in range(1, max_dist + 1):
        i = np.flatnonzero(rest >= d)
        key = 4 * site[i] + 2 * meth[i] + meth[i + d]
        tables[d] = np.bincount(key, weights=weight[i], minlength=4 * n_sites).reshape(n_sites, 4)
    return tables


def unlinked(table, r2_cut, p_cut):
    """Whether CpG pairs fail the linkage test, given their 00/01/10/11 counts.

    Same arithmetic as Util.getR2FromMap: r2 is NaN when either site is monomorphic, and
    a NaN r2 never falls below the cutoff (such pairs fail on the p-value instead).
    """
    n00, n01, n10, n11 = table.T
    n = n00 + n01 + n10 + n11
    covered = n > 0
    with np.errstate(divide='ignore', invalid='ignore'):
        pa = (n10 + n11) / n
        pb = (n01 + n11) / n
        d = n11 / n - pa * pb
        den = pa * (1 - pa) * pb * (1 - pb)
        r2 = np.where(den == 0, np.nan, d * d / den)
        r2 = np.where(d < 0, -r2, r2)
        pval = np.ones_like(n)
        pval[covered] = binom.sf(n11[covered] - 1, n[covered], (pa * pb)[covered])
        return ~covered | (r2 < r2_cut) | (pval > p_cut)


def scan_blocks(fails, n, window):
    """Blocks (first, last) of site indices, following the scan in MHBDiscovery.getMHB.

    fails[d][k] tells whether the pair (k, k + d) is unlinked, for k = 0 .. n - d - 1.
    """
    # distance to the nearest site before e that e is not linked with (0: linked to all in window)
    break_dist = np.zeros(n, dtype=np.int64)
    for d in range(window - 1, 0, -1):
        if d < n:
            break_dist[d:][fails[d]] = d
    ends = np.flatnonzero(break_dist)
    if ends.size == 0:
        return [(0, n - 1)] if n - 1 >= window else []
    new_start = np.maximum.accumulate(ends - break_dist[ends] + 1)
    start = np.concatenate(([0], new_start[:-1]))  # start of the block each break closes
    keep = ends - start >= window
    blocks = list(zip(start[keep].tolist(), (ends[keep] - 1).tolist()))
    if n - 1 - new_start[-1] >= window:
        blocks.append((int(new_start[-1]), n - 1))
    return blocks


def chunk_mhbs(mhap, cpg, chrom, start, end, window, r2_cut, p_cut, stats):
    """MHBs (start, end) found in one chunk [start, end]."""
    starts, ends, haps, counts = fetch_mhaps(mhap, chrom, start, end)
    if not haps:
        return []
    # CpG sites spanned by the mHaps, to locate each call; blocks use only sites inside the chunk
    pos = fetch_cpgs(cpg, chrom, min(start, starts.min()), max(end, ends.max()))
    lo, hi = np.searchsorted(pos, start, 'left'), np.searchsorted(pos, end, 'right')
    n = hi - lo
    if n < 1:
        return []

    # each mHap must start on a CpG site and span len(hap) consecutive sites
    lens = np.fromiter(map(len, haps), dtype=np.int64, count=len(haps))
    first = np.searchsorted(pos, starts)
    last = first + lens - 1
    ok = (last < len(pos)) & (pos[np.minimum(first, len(pos) - 1)] == starts) & \
         (pos[np.minimum(last, len(pos) - 1)] == ends)
    if not ok.all():
        stats['skipped'] += int((~ok).sum())
        haps = [h for h, keep in zip(haps, ok) if keep]
        first, lens, counts = first[ok], lens[ok], counts[ok]

    max_dist = min(window - 1, n - 1)
    tables = pair_counts(first, haps, lens, counts, len(pos), max_dist)
    fails = {d: unlinked(tables[d][lo:hi - d], r2_cut, p_cut) for d in range(1, max_dist + 1)}
    sites = pos[lo:hi]
    return [(int(sites[a]), int(sites[b])) for a, b in scan_blocks(fails, n, window)]


def discover_mhbs(mhap_path, cpg_path, regions, window=5, r2_cut=0.5, p_cut=0.05, stats=None):
    """Yield the MHBs (chrom, start, end) found in regions given as 1-based closed intervals."""
    if stats is None:
        stats = {}
    stats.setdefault('skipped', 0)
    mhap, cpg = TabixFile(mhap_path), TabixFile(cpg_path)
    seen = set()
    for chrom, start, end in regions:
        for chunk_start, chunk_end in split_region(start, end):
            for mhb_start, mhb_end in chunk_mhbs(mhap, cpg, chrom, chunk_start, chunk_end,
                                                 window, r2_cut, p_cut, stats):
                key = (chrom, mhb_start, mhb_end)
                if key not in seen:
                    seen.add(key)
                    yield key


def mhb_discovery(mhap_path, cpg_path, output_file, region=None, bed_path=None,
                  window=5, r_square=0.5, p_value=0.05):
    """Write the MHBs of a region, of the regions in a BED file, or of the whole genome."""
    if region is not None and bed_path is not None:
        raise ValueError('--region and --bedPath can not be given at the same time')
    if window < 2:
        raise ValueError('--window must be at least 2')
    out_dir = os.path.dirname(output_file)
    if out_dir and not os.path.isdir(out_dir):
        raise FileNotFoundError('output directory does not exist: ' + out_dir)

    t0 = time.time()
    if region is not None:
        regions = [parse_region(region)]
    elif bed_path is not None:
        regions = read_bed(bed_path)
    else:
        regions = genome_regions(cpg_path)

    stats = {}
    n_mhb = 0
    with open(output_file, 'w') as out:
        for chrom, start, end in discover_mhbs(mhap_path, cpg_path, regions, window, r_square, p_value, stats):
            out.write('%s\t%d\t%d\n' % (chrom, start, end))
            n_mhb += 1
    if stats['skipped']:
        print('Warning: %d mHaps skipped because their CpG positions do not match %s'
              % (stats['skipped'], cpg_path), file=sys.stderr)
    print('%d MHBs written to %s, MHB time span: %.1f s' % (n_mhb, output_file, time.time() - t0))
    return n_mhb


def add_arguments(parser):
    parser.add_argument('--mhapPath', type=str, required=True, help='input file, mhap.gz format, generated by mHapTools and indexed')
    parser.add_argument('--cpgPath', type=str, required=True, help='genomic CpG file, gz format and indexed')
    where = parser.add_mutually_exclusive_group()
    where.add_argument('--region', type=str, default=None, help='one region, in the format of chr:start-end')
    where.add_argument('--bedPath', type=str, default=None, help='input BED file; without --region and --bedPath, the whole genome is scanned')
    parser.add_argument('--outputFile', type=str, required=True, help='output file name')
    parser.add_argument('--window', type=int, default=5, help='size of core window [5]')
    parser.add_argument('--r_square', type=float, default=0.5, help='R-square cutoff [0.5]')
    parser.add_argument('--p_value', type=float, default=0.05, help='P-value cutoff [0.05]')


def run(args):
    mhb_discovery(args.mhapPath, args.cpgPath, args.outputFile, region=args.region, bed_path=args.bedPath,
                  window=args.window, r_square=args.r_square, p_value=args.p_value)


def main(argv=None):
    parser = argparse.ArgumentParser(description='identification of methylation haplotype blocks within a region or genome-wide')
    add_arguments(parser)
    run(parser.parse_args(argv))


if __name__ == '__main__':
    main()

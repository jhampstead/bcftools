#!/usr/bin/env python3
#
#   The GERP scores for bcftools +lof --gerp, extracted from a bigWig file such as the GERP file of
#   LOFTEE, gerp_conservation_scores.homo_sapiens.GRCh38.bw:
#
#       lof-gerp.py -g genes.gff3.gz gerp.bw | bgzip > gerp.tsv.gz
#       tabix -p bed gerp.tsv.gz
#       bcftools +lof in.vcf.gz -- -f ref.fa -g genes.gff3.gz ... --gerp gerp.tsv.gz
#
#   LOFTEE reads the scores through the bigWig library of Jim Kent, which averages them from the
#   file's zoom levels for all but the shortest intervals, so its GERP-weighted distances depend on
#   how the bigWig file was written. To give the same distances, +lof needs the zoom records as well
#   as the scores, for the exons it can query: those of the transcripts with a CDS in the GFF3 file
#   given to bcftools csq and +lof. Only the Python standard library is needed.
#
#   The output is tab-delimited with the columns chrom, beg, end (0-based, half-open), level,
#   count and value, sorted for tabix -p bed:
#       level -1    a region covered by the file: an exon, or exons merged
#       level 0     scores of the bigWig file: count bases each from beg, values separated by ','
#       level N     a zoom record with reduction level N: count bases with data, value their sum
#   The header lists the zoom levels and the chromosomes of the bigWig file.
#
import argparse
import array
import bisect
import gzip
import os
import struct
import sys
import zlib

BIGWIG_MAGIC = 0x888FFC26
BPT_MAGIC = 0x78CA8C91
CIRTREE_MAGIC = 0x2468ACE0
BEDGRAPH, VARSTEP, FIXEDSTEP = 1, 2, 3


class BigWig:
    """The parts of a bigWig file that LOFTEE's queries use, after the format of Jim Kent's library
    (bbiFile.h, bwgQuery.c, bPlusTree.c, cirTree.c)."""

    def __init__(self, path):
        self.fh = open(path, 'rb')
        head = self.fh.read(64)
        if len(head) == 64 and struct.unpack('<I', head[:4])[0] == BIGWIG_MAGIC:
            self.e = '<'
        elif len(head) == 64 and struct.unpack('>I', head[:4])[0] == BIGWIG_MAGIC:
            self.e = '>'
        else:
            sys.exit(f'Error: {path} is not a bigWig file')
        (_, _, nzoom, chrom_tree, _, self.full_index, _, _, _, _, self.uncompress, _) = \
            struct.unpack(self.e + 'IHHQQQHHQQIQ', head)
        self.zooms = []         # (reduction level, R-tree offset), in file order
        data = self.fh.read(24 * nzoom)
        for i in range(nzoom):
            reduction, _, _, index = struct.unpack_from(self.e + 'IIQQ', data, 24 * i)
            self.zooms.append((reduction, index))
        self.chroms = {}        # name -> (id, size), in file order
        self.read_chroms(chrom_tree)

    def read(self, offset, size):
        self.fh.seek(offset)
        return self.fh.read(size)

    def read_chroms(self, offset):
        magic, _, key_size, val_size, _, _ = struct.unpack(self.e + 'IIIIQQ', self.read(offset, 32))
        if magic != BPT_MAGIC: sys.exit('Error: the chromosome tree of the bigWig file is damaged')
        nodes = [offset + 32]
        while nodes:
            node = nodes.pop(0)
            leaf, _, count = struct.unpack(self.e + 'BBH', self.read(node, 4))
            item = key_size + (val_size if leaf else 8)
            data = self.fh.read(item * count)
            for i in range(count):
                key = data[i * item:i * item + key_size].rstrip(b'\0').decode()
                if leaf:
                    self.chroms[key] = struct.unpack_from(self.e + 'II', data, i * item + key_size)
                else:
                    nodes.append(struct.unpack_from(self.e + 'Q', data, i * item + key_size)[0])

    def blocks(self, index, chrom_id):
        """The data blocks (offset, size) of the R-tree at index that can hold data of a chromosome."""
        magic = struct.unpack(self.e + 'I', self.read(index, 4))[0]
        if magic != CIRTREE_MAGIC: sys.exit('Error: an index of the bigWig file is damaged')
        out, nodes = [], [index + 48]
        while nodes:
            node = nodes.pop(0)
            leaf, _, count = struct.unpack(self.e + 'BBH', self.read(node, 4))
            if leaf:
                for sc, _, ec, _, offset, size in struct.iter_unpack(self.e + 'IIIIQQ', self.fh.read(32 * count)):
                    if sc <= chrom_id <= ec: out.append((offset, size))
            else:
                for sc, _, ec, _, offset in struct.iter_unpack(self.e + 'IIIIQ', self.fh.read(24 * count)):
                    if sc <= chrom_id <= ec: nodes.append(offset)
        return sorted(out)

    def block(self, offset, size):
        data = self.read(offset, size)
        return zlib.decompress(data) if self.uncompress else data

    def numbers(self, typecode, data):
        a = array.array(typecode, data)
        if (self.e == '<') != (sys.byteorder == 'little'): a.byteswap()
        return a

    def section(self, data):
        """The items of a data block, one section: chromosome id and the arrays of item starts, ends
        and values."""
        chrom_id, start, _, step, span, kind, _, count = struct.unpack_from(self.e + 'IIIIIBBH', data, 0)
        data = data[24:]
        if kind == BEDGRAPH:
            ints, vals = self.numbers('I', data), self.numbers('f', data)
            return chrom_id, ints[0::3], ints[1::3], vals[2::3]
        if kind == VARSTEP:
            ints, vals = self.numbers('I', data), self.numbers('f', data)
            starts = ints[0::2]
            return chrom_id, starts, array.array('I', (s + span for s in starts)), vals[1::2]
        if kind == FIXEDSTEP:
            starts = array.array('I', range(start, start + step * count, step)) if count else array.array('I')
            return chrom_id, starts, array.array('I', (s + span for s in starts)), self.numbers('f', data[:4 * count])
        sys.exit(f'Error: unknown bigWig section type {kind}')

    def zoom_records(self, data):
        """The zoom records of a block: chromosome ids, starts, ends, valid counts and sums."""
        ints, vals = self.numbers('I', data), self.numbers('f', data)
        return ints[0::8], ints[1::8], ints[2::8], ints[3::8], vals[6::8]


_float32 = {}
def fmt(x):
    """The shortest decimal that reads back as the same 32-bit float."""
    s = _float32.get(x)
    if s is None:
        for p in range(1, 10):
            s = '%.*g' % (p, x)
            if struct.unpack('f', struct.pack('f', float(s)))[0] == x: break
        _float32[x] = s
    return s


def find_chr(name, names):
    """The name of a chromosome in names, as +lof finds it: as is, then with or without 'chr', and
    MT for chrM and back."""
    if name in names: return name
    alt = 'chrM' if name == 'MT' else 'MT' if name == 'chrM' else name[3:] if name.startswith('chr') else 'chr' + name
    return alt if alt in names else None


def read_exons(path):
    """The exons of the transcripts with a CDS, by chromosome: 0-based, half-open, merged."""
    exons, coding = {}, set()
    with open(path, 'rb') as fh:
        gz = fh.read(2) == b'\x1f\x8b'
    with (gzip.open(path, 'rt') if gz else open(path)) as fh:
        for line in fh:
            if line.startswith('#'): continue
            col = line.rstrip('\n').split('\t')
            if len(col) < 9 or col[2] not in ('exon', 'CDS'): continue
            parents = next((a[7:] for a in col[8].split(';') if a.startswith('Parent=')), None)
            if not parents: continue
            for parent in parents.split(','):
                if col[2] == 'CDS': coding.add(parent)
                else: exons.setdefault(parent, []).append((col[0], int(col[3]) - 1, int(col[4])))
    regions = {}
    for tscript in coding:
        for chrom, beg, end in exons.get(tscript, []):
            regions.setdefault(chrom, []).append((beg, end))
    for chrom, reg in regions.items():
        reg.sort()
        merged = [list(reg[0])]
        for beg, end in reg[1:]:
            if beg <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], end)
            else: merged.append([beg, end])
        regions[chrom] = [tuple(r) for r in merged]
    return regions


def overlapping(starts, ends, reg_beg, reg_end):
    """Index ranges of the sorted, non-overlapping items overlapping any of the regions, merged."""
    ranges = []
    if not len(starts): return ranges
    first = bisect.bisect_right(reg_end, starts[0])
    for r in range(first, len(reg_beg)):
        if reg_beg[r] >= ends[-1]: break
        i = bisect.bisect_right(ends, reg_beg[r])          # the first item ending after the region start
        j = bisect.bisect_left(starts, reg_end[r])         # the first item starting at or after its end
        if i >= j: continue
        if ranges and i <= ranges[-1][1]: ranges[-1][1] = max(ranges[-1][1], j)
        else: ranges.append([i, j])
    return ranges


def export_chrom(bw, name, regions, out):
    chrom_id = bw.chroms[name][0]
    reg_beg, reg_end = [r[0] for r in regions], [r[1] for r in regions]
    lines = [(beg, f'{name}\t{beg}\t{end}\t-1\t.\t.') for beg, end in regions]
    nscores = 0

    for offset, size in bw.blocks(bw.full_index, chrom_id):
        cid, starts, ends, vals = bw.section(bw.block(offset, size))
        if cid != chrom_id: continue
        for i, j in overlapping(starts, ends, reg_beg, reg_end):
            # runs of contiguous items of the same length, one line each
            k = i
            while k < j:
                span, m = ends[k] - starts[k], k + 1
                while m < j and starts[m] == ends[m - 1] and ends[m] - starts[m] == span: m += 1
                lines.append((starts[k], f'{name}\t{starts[k]}\t{ends[m - 1]}\t0\t{span}\t' +
                              ','.join(fmt(v) for v in vals[k:m])))
                nscores += m - k
                k = m

    nzoom = 0
    for reduction, index in bw.zooms:
        for offset, size in bw.blocks(index, chrom_id):
            cids, starts, ends, valid, sums = bw.zoom_records(bw.block(offset, size))
            # a block can hold records of more than one chromosome
            keep = [k for k in range(len(cids)) if cids[k] == chrom_id]
            if not keep: continue
            if len(keep) < len(cids):
                starts, ends = array.array('I', (starts[k] for k in keep)), array.array('I', (ends[k] for k in keep))
                valid, sums = [valid[k] for k in keep], [sums[k] for k in keep]
            for i, j in overlapping(starts, ends, reg_beg, reg_end):
                for k in range(i, j):
                    lines.append((starts[k], f'{name}\t{starts[k]}\t{ends[k]}\t{reduction}\t{valid[k]}\t{fmt(sums[k])}'))
                nzoom += j - i

    lines.sort(key=lambda x: x[0])
    for _, line in lines:
        out.write(line + '\n')
    return nscores, nzoom


def main():
    parser = argparse.ArgumentParser(
        description='Extract from a GERP bigWig file the scores and zoom records that bcftools +lof --gerp '
                    'needs to compute the GERP-weighted distances of LOFTEE exactly. Writes to standard output, '
                    'compress it with bgzip and index it with "tabix -p bed".',
        epilog='Example: lof-gerp.py -g genes.gff3.gz gerp.bw | bgzip > gerp.tsv.gz && tabix -p bed gerp.tsv.gz')
    parser.add_argument('-g', '--gff', required=True, help='the GFF3 file given to bcftools csq and +lof')
    parser.add_argument('bigwig', help='the GERP bigWig file, e.g. gerp_conservation_scores.homo_sapiens.GRCh38.bw')
    args = parser.parse_args()

    bw = BigWig(args.bigwig)
    regions = read_exons(args.gff)
    out = sys.stdout
    out.write('##lof_gerp=1\n')
    out.write(f'##bigwig={os.path.basename(args.bigwig)}\n')
    out.write(f'##gff={os.path.basename(args.gff)}\n')
    out.write('##zoom_levels=' + ','.join(str(z[0]) for z in bw.zooms) + '\n')
    for name, (_, size) in bw.chroms.items():
        out.write(f'##chrom={name}:{size}\n')
    out.write('#chrom\tbeg\tend\tlevel\tcount\tvalue\n')

    nscores = nzoom = 0
    for gff_chrom in sorted(regions):
        name = find_chr(gff_chrom, bw.chroms)
        if name is None:
            print(f'Warning: {gff_chrom} is not in {args.bigwig}, skipped {len(regions[gff_chrom])} regions',
                  file=sys.stderr)
            continue
        n, z = export_chrom(bw, name, regions[gff_chrom], out)
        nscores += n
        nzoom += z
    print(f'Exported {nscores} scores and {nzoom} zoom records for {sum(len(r) for r in regions.values())} '
          f'regions on {len(regions)} chromosomes', file=sys.stderr)


if __name__ == '__main__':
    main()

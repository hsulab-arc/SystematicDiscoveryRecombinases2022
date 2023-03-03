from umi_tools import UMIClusterer
from pybedtools import BedTool
from Bio import pairwise2
import pysam
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
import gzip
from snakemake import shell
import os


def fastq_gen(fq):
    with gzip.open(fq, 'rt') as handle:
        while True:
            l1 = handle.readline().strip()
            l2 = handle.readline().strip()
            l3 = handle.readline().strip()
            l4 = handle.readline().strip()
            if len(l2) == 0:
                break
            yield [l1, l2, l3, l4]

def revcomp(seq): # KEEP
    seq = seq.upper()
    out = ''
    for c in seq[::-1]:
        if c == 'A':
            out += 'T'
        elif c == 'T':
            out += 'A'
        elif c == 'C':
            out += 'G'
        elif c == 'G':
            out += 'C'
        else:
            out += c
    return out

def get_kmers(s, k=8):
    s = s.upper()
    for i in range(len(s) - k + 1):
        yield s[i:i + k]

def clean_umi(input, output, params): # KEEP

    metadata = pd.read_csv(params.metadata, sep='\t')
    stagger_seq = list(metadata[metadata['sample'] == params.sample]['i7_stagger'])[0]
    try:
        stagger_length = len(stagger_seq)
    except:
        stagger_length = 0

    primer_seq = list(metadata[metadata['sample'] == params.sample]['primer_seq'])[0]
    umi_start = int(list(metadata[metadata['sample'] == params.sample]['umi_start'])[0])
    umi_end = int(list(metadata[metadata['sample'] == params.sample]['umi_end'])[0])

    print(primer_seq)
    print(umi_start, umi_end)

    fq1 = fastq_gen(input.r1)
    fq2 = fastq_gen(input.r2)

    outfq1 = gzip.open(output.r1, 'wt')
    outfq2 = gzip.open(output.r2, 'wt')

    for rec1, rec2 in zip(fq1, fq2):

        primer_start = rec2[1].find(primer_seq)
        if primer_start != -1:
            rec2[1] = rec2[1][primer_start:]
            rec2[3] = rec2[3][primer_start:]
        else:
            rec2[1] = rec2[1][stagger_length:]
            rec2[3] = rec2[3][stagger_length:]

        umi_seq = rec2[1][umi_start:umi_end]
        umi_seq_revcomp = revcomp(umi_seq)
        s1, s2, s3 = rec2[1][:umi_start], umi_seq, rec2[1][umi_end:]
        s2 = 'A' * len(umi_seq)

        rec2[1] = s1 + s2 + s3
        rec1[1] = rec1[1].replace(umi_seq_revcomp, 'T' * len(umi_seq))

        rec1[0] = rec1[0].split()
        rec1[0][0] = rec1[0][0] + ':UMI_' + umi_seq
        rec1[0] = ' '.join(rec1[0])

        rec2[0] = rec2[0].split()
        rec2[0][0] = rec2[0][0] + ':UMI_' + umi_seq
        rec2[0] = ' '.join(rec2[0])

        print(*rec1, sep='\n', file=outfq1)
        print(*rec2, sep='\n', file=outfq2)

    outfq1.close()
    outfq2.close()

def get_perc_identity(read): # KEEP
    total = read.query_alignment_length
    identity = (total - read.get_tag('NM')) / total
    return identity

def add_read_info(bam, label, read_info): # KEEP
    print("Getting read info for", label)
    for read in bam:

        if read.is_supplementary:
            continue
        if read.is_unmapped:
            read_info[read.query_name]
            continue

        refstart, refend = read.reference_start, read.reference_end
        orient = read.is_reverse
        if orient == True:
            orient = '-'
        else:
            orient = '+'

        clipleft, clipright = 0, 0
        cigtups = read.cigartuples
        if cigtups[0][0] == 4:
            clipleft = cigtups[0][1]
        if cigtups[-1][0] == 4:
            clipright = cigtups[-1][1]
        out = {'refname': read.reference_name, 'refstart': refstart, 'refend': refend, 'orient': orient,
               'clipleft': clipleft, 'clipright': clipright, 'pident': get_perc_identity(read)}
        read_info[read.query_name][label] = out

def format_outline(outfile, read, r1_category, r2_category, primer_status, r1_attd_status, r2_attd_status,
                   r1_human_info, r2_human_info, r1_donor_info, r2_donor_info, r1_genome_blacklist,
                   r2_genome_blacklist, r1_donor_check, r2_donor_check): # KEEP
    if r1_human_info is not None:
        r1_human_info = r1_human_info['refname'], r1_human_info['refstart'], r1_human_info['refend'], r1_human_info[
            'orient'], r1_human_info['clipleft'], r1_human_info['clipright'], r1_human_info[
                            'pident'], r1_genome_blacklist
    else:
        r1_human_info = ["NA"] * 8

    if r2_human_info is not None:
        r2_human_info = r2_human_info['refname'], r2_human_info['refstart'], r2_human_info['refend'], r2_human_info[
            'orient'], r2_human_info['clipleft'], r2_human_info['clipright'], r2_human_info[
                            'pident'], r2_genome_blacklist
    else:
        r2_human_info = ["NA"] * 8

    if r1_donor_info is not None:
        r1_donor_info = r1_donor_info['refname'], r1_donor_info['refstart'], r1_donor_info['refend'], r1_donor_info[
            'orient'], r1_donor_info['clipleft'], r1_donor_info['clipright'], r1_donor_info['pident'],
    else:
        r1_donor_info = ["NA"] * 7

    if r2_donor_info is not None:
        r2_donor_info = r2_donor_info['refname'], r2_donor_info['refstart'], r2_donor_info['refend'], r2_donor_info[
            'orient'], r2_donor_info['clipleft'], r2_donor_info['clipright'], r2_donor_info['pident'],
    else:
        r2_donor_info = ["NA"] * 7

    print(
        read, r1_category, r2_category, primer_status, r1_attd_status, r2_attd_status, *r1_human_info, *r2_human_info,
        *r1_donor_info, *r2_donor_info, r1_donor_check, r2_donor_check, sep='\t', file=outfile
    )

def get_genome_blacklist(human_info, blacklist): # KEEP
    out_blacklist = 'None'
    if human_info['refname'] in blacklist:
        for r in blacklist[human_info['refname']]:
            if human_info['refstart'] >= r['start'] and human_info['refstart'] <= r['end']:
                out_blacklist = r['name']
    return out_blacklist

def run_donor_check(bam, correct_donor):
    top_donor = defaultdict(lambda : {'best_score': None, 'top_donors': None})
    for read in bam:
        pident, alnlen = get_perc_identity(read), read.reference_end - read.reference_start
        score = pident*alnlen
        if top_donor[read.query_name]['best_score'] is None or top_donor[read.query_name]['best_score'] < score:
            top_donor[read.query_name]['best_score'] = score
            top_donor[read.query_name]['top_donors'] = set()
            for refname in read.reference_name.split('|'):
                top_donor[read.query_name]['top_donors'].add(refname)
        elif top_donor[read.query_name]['best_score'] == score:
            for refname in read.reference_name.split('|'):
                top_donor[read.query_name]['top_donors'].add(refname)

    has_correct_donor = dict()
    for read in top_donor:
        if correct_donor in top_donor[read]['top_donors']:
            has_correct_donor[read] = True
        else:
            has_correct_donor[read] = False

    return has_correct_donor


def analyze_reads(input, output, params): # KEEP

    print("Collecting all read names..")
    all_reads = set()
    for l1, l2, l3, l4 in fastq_gen(input.r1):
        all_reads.add(l1.split()[0][1:])

    # Collecting donor check information
    print("Getting donor check info for R1")
    donor_check_r1 = run_donor_check(pysam.AlignmentFile(input.donor_check_bam1, 'rb'), params.sample)
    print("Getting donor check info for R2")
    donor_check_r2 = run_donor_check(pysam.AlignmentFile(input.donor_check_bam2, 'rb'), params.sample)

    read_info = defaultdict(lambda: {"R1_human": None, "R2_human": None, "R1_donor": None, "R2_donor": None})
    # Collecting read information
    human_bam1 = pysam.AlignmentFile(input.human_bam1, 'rb')
    add_read_info(human_bam1, 'R1_human', read_info)

    human_bam2 = pysam.AlignmentFile(input.human_bam2, 'rb')
    add_read_info(human_bam2, 'R2_human', read_info)

    donor_bam1 = pysam.AlignmentFile(input.donor_bam1, 'rb')
    add_read_info(donor_bam1, 'R1_donor', read_info)

    donor_bam2 = pysam.AlignmentFile(input.donor_bam2, 'rb')
    add_read_info(donor_bam2, 'R2_donor', read_info)

    # Getting primer and attd coordinates
    attd_start, attd_end = None, None
    primer_start = None
    with open(input.bed) as infile:
        for line in infile:
            line = line.strip().split()
            if line[3] == "Primer":
                primer_start = int(line[1])
            elif line[3] == 'attD':
                attd_start, attd_end = int(line[1]), int(line[2])

    # Getting blacklisted regions
    blacklist = defaultdict(list)
    with open(params.blacklist) as infile:
        for line in infile:
            line = line.strip().split()
            blacklist[line[0]].append({'start': int(line[1]), 'end': int(line[2]), 'name': line[3]})

    # Writing results to file
    outfile = gzip.open(output.summary, 'wt')
    print("query_name", "r1_align", "r2_align", "primer_status", "r1_attd_status", "r2_attd_status", 'r1_human_refname',
          'r1_human_refstart', 'r1_human_refend', 'r1_human_orient', 'r1_human_clipleft', 'r1_human_clipright',
          'r1_human_pident', 'r1_genome_blacklist', 'r2_human_refname', 'r2_human_refstart', 'r2_human_refend',
          'r2_human_orient', 'r2_human_clipleft', 'r2_human_clipright', 'r2_human_pident', 'r2_genome_blacklist',
          'r1_donor_refname', 'r1_donor_refstart', 'r1_donor_refend', 'r1_donor_orient', 'r1_donor_clipleft',
          'r1_donor_clipright', 'r1_donor_pident', 'r2_donor_refname', 'r2_donor_refstart', 'r2_donor_refend',
          'r2_donor_orient', 'r2_donor_clipleft', 'r2_donor_clipright', 'r2_donor_pident', 'r1_donor_check',
          'r2_donor_check', sep='\t', file=outfile)

    for read in all_reads:

        info = read_info[read]
        r1_donor_info, r2_donor_info = info['R1_donor'], info['R2_donor']
        r1_human_info, r2_human_info = info['R1_human'], info['R2_human']

        maps_r1_human = r1_human_info is not None
        maps_r2_human = r2_human_info is not None
        maps_r1_donor = r1_donor_info is not None
        maps_r2_donor = r2_donor_info is not None

        r1_genome_blacklist = 'NA'
        r2_genome_blacklist = 'NA'

        r1_category = 'NOTHING'
        if maps_r1_human and maps_r1_donor:
            r1_category = 'GENOMIC+DONOR'
        elif maps_r1_human:
            r1_category = 'GENOMIC_ONLY'
        elif maps_r1_donor:
            r1_category = 'DONOR_ONLY'

        if 'GENOMIC' in r1_category:
            r1_genome_blacklist = get_genome_blacklist(r1_human_info, blacklist)

        r2_category = 'NOTHING'
        if maps_r2_human and maps_r2_donor:
            r2_category = 'GENOMIC+DONOR'
        elif maps_r2_human:
            r2_category = 'GENOMIC_ONLY'
        elif maps_r2_donor:
            r2_category = 'DONOR_ONLY'

        if 'GENOMIC' in r2_category:
            r2_genome_blacklist = get_genome_blacklist(r2_human_info, blacklist)

        primer_status = 'NA'
        if maps_r2_donor:
            if primer_start - 2 <= r2_donor_info['refstart'] <= primer_start + 2:
                primer_status = 'PRIMER'
            else:
                primer_status = 'NOPRIMER'

        r1_attd_status = 'NA'
        if maps_r1_donor:
            if r1_donor_info['refstart'] < attd_start and r1_donor_info['refend'] > attd_end:
                r1_attd_status = 'INTACT'
            elif attd_start < r1_donor_info['refstart'] < attd_end < r1_donor_info['refend'] and r1_donor_info['clipleft'] > 4:
                r1_attd_status = 'CLIPPED_LEFT'
            elif r1_donor_info['refend'] > attd_start > r1_donor_info['refstart'] < attd_end and r1_donor_info['clipright'] > 4:
                r1_attd_status = 'CLIPPED_RIGHT'
            else:
                r1_attd_status = 'UNKNOWN'

        r2_attd_status = 'NA'
        if maps_r2_donor:
            if r2_donor_info['refstart'] < attd_start and r2_donor_info['refend'] > attd_end:
                r2_attd_status = 'INTACT'
            elif attd_start < r2_donor_info['refstart'] < attd_end < r2_donor_info['refend'] and r2_donor_info['clipleft'] > 4:
                r2_attd_status = 'CLIPPED_LEFT'
            elif r2_donor_info['refend'] > attd_start > r2_donor_info['refstart'] < attd_end and r2_donor_info['clipright'] > 4:
                r2_attd_status = 'CLIPPED_RIGHT'
            else:
                r2_attd_status = 'UNKNOWN'

        if read not in donor_check_r1:
            r1_donor_check = 'UNKNOWN'
        else:
            r1_donor_check = donor_check_r1[read]
        if read not in donor_check_r2:
            r2_donor_check = 'UNKNOWN'
        else:
            r2_donor_check = donor_check_r2[read]

        format_outline(outfile, read, r1_category, r2_category, primer_status, r1_attd_status, r2_attd_status,
                       r1_human_info, r2_human_info, r1_donor_info, r2_donor_info, r1_genome_blacklist,
                       r2_genome_blacklist, r1_donor_check, r2_donor_check)

    outfile.close()

def get_junction_reads(input, output): # KEEP

    junction_reads = set()
    with gzip.open(input.summary, 'rt') as infile:
        header = infile.readline().strip().split('\t')
        for line in infile:
            line = line.strip().split('\t')
            line = {header[i]: line[i] for i in range(len(line))}
            if line['r1_donor_check'] == 'False' or line['r2_donor_check'] != 'True':
                continue
            if (line['r1_align'] == 'GENOMIC_ONLY' and line['r2_align'] == 'GENOMIC_ONLY'):
                continue
            if (line['r1_align'] == 'DONOR_ONLY' and line['r2_align'] == 'DONOR_ONLY'):
                continue
            if ("ONLY" in line['r1_align'] and line['r2_align'] == 'NOTHING') or (
                    "ONLY" in line['r2_align'] and line['r1_align'] == 'NOTHING'):
                continue
            if line['r2_attd_status'] == 'INTACT' or line['r1_attd_status'] == 'INTACT':
                continue
            if line['primer_status'] != 'PRIMER':
                continue
            if line['r1_genome_blacklist'] != 'None' and line['r1_genome_blacklist'] != 'NA':
                continue
            if line['r2_genome_blacklist'] != 'None' and line['r2_genome_blacklist'] != 'NA':
                continue
            if 'GENOMIC' not in line['r1_align']:
                continue
            junction_reads.add(line['query_name'])

    out1 = gzip.open(output.r1, 'wt')
    out2 = gzip.open(output.r2, 'wt')

    for r1, r2 in zip(fastq_gen(input.r1), fastq_gen(input.r2)):
        query_name = r1[0].split()[0].replace('@', '')
        if query_name not in junction_reads:
            continue
        print(*r1, sep='\n', file=out1)
        print(*r2, sep='\n', file=out2)

    out1.close()
    out2.close()

def get_donor_reads(input, output):
    donor_reads = set()
    with gzip.open(input.summary, 'rt') as infile:
        header = infile.readline().strip().split('\t')
        for line in infile:
            line = line.strip().split()
            line = {header[i]: line[i] for i in range(len(line))}
            if line['r1_donor_check'] == 'False' or line['r2_donor_check'] != 'True':
                continue
            if not (line['r1_align'] == 'DONOR_ONLY' and line['r2_align'] == 'DONOR_ONLY'):
                continue
            if line['r1_attd_status'] == 'CLIPPED_LEFT' or line['r2_attd_status'] == 'CLIPPED_LEFT':
                continue
            if line['r1_attd_status'] == 'CLIPPED_RIGHT' or line['r2_attd_status'] == 'CLIPPED_RIGHT':
                continue
            if line['primer_status'] != 'PRIMER':
                continue
            donor_reads.add(line['query_name'])

    out1 = gzip.open(output.r1, 'wt')
    out2 = gzip.open(output.r2, 'wt')

    for r1, r2 in zip(fastq_gen(input.r1), fastq_gen(input.r2)):
        query_name = r1[0].split()[0].replace('@', '')
        if query_name not in donor_reads:
            continue
        print(*r1, sep='\n', file=out1)
        print(*r2, sep='\n', file=out2)

    out1.close()
    out2.close()

def get_percent_ident(seq1, seq2):
    orig_seq1, orig_seq2 = seq1, seq2
    orig_nogaps_seq1, orig_nogaps_seq2 = orig_seq1.replace('-', ''), orig_seq2.replace('-', '')

    trim_left = max([len(seq1) - len(seq1.lstrip('-')), len(seq2) - len(seq2.lstrip('-'))])
    seq1, seq2 = seq1[trim_left:], seq2[trim_left:]

    trim_right = max([len(seq1) - len(seq1.rstrip('-')), len(seq2) - len(seq2.rstrip('-'))])
    if trim_right != 0:
        seq1, seq2 = seq1[:-trim_right], seq2[:-trim_right]

    seq1_tmp, seq2_tmp = '', ''
    for i in range(len(seq1)):
        if not (seq1[i] == '-' and seq2[i] == '-'):
            seq1_tmp += seq1[i]
            seq2_tmp += seq2[i]
    seq1, seq2 = seq1_tmp, seq2_tmp

    matches = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            matches += 1

    seq1_gapped = ''
    seq2_gapped = ''
    for i in range(len(seq1)):

        if seq1[i] == '-' or seq2[i] == '-' or seq1[i] == 'N' or seq2[i] == 'N':
            seq1_gapped += '-'
            seq2_gapped += '-'
        else:
            seq1_gapped += seq1[i]
            seq2_gapped += seq2[i]

    seq1_gapped = seq1_gapped.split('-')
    seq2_gapped = seq2_gapped.split('-')

    longest_stretch_seq1 = ''
    longest_stretch_seq2 = ''
    contig_matches = 0
    longest_stretch = None
    for c1, c2 in zip(seq1_gapped, seq2_gapped):

        if longest_stretch is None:
            longest_stretch_seq1, longest_stretch_seq2 = c1, c2
            longest_stretch = len(c1)
            this_contig_matches = 0
            for i in range(len(c1)):
                if c1[i] == c2[i]:
                    this_contig_matches += 1
            contig_matches = this_contig_matches

        elif len(c1) > longest_stretch:
            longest_stretch_seq1, longest_stretch_seq2 = c1, c2
            longest_stretch = len(c1)
            this_contig_matches = 0
            for i in range(len(c1)):
                if c1[i] == c2[i]:
                    this_contig_matches += 1
            contig_matches = this_contig_matches

    return matches / len(seq1), len(seq1), contig_matches / longest_stretch, longest_stretch, orig_nogaps_seq2.find(
        seq2), orig_nogaps_seq2.find(seq2) + len(seq2), orig_nogaps_seq2.find(
        longest_stretch_seq2), orig_nogaps_seq2.find(longest_stretch_seq2) + len(longest_stretch_seq2)

def get_align_core(seq, attd, pos, orient, core_start, core_end): # KEEP
    nogaps_seq = seq.replace('-', '')

    aln_core_start, aln_core_end = None, None
    attd_pos = -1
    for i in range(len(attd)):
        if attd[i] == '-':
            continue
        attd_pos += 1
        if attd_pos == core_start:
            aln_core_start = i
        elif attd_pos == core_end:
            aln_core_end = i

    # print(seq)
    # print(attd)
    genomic_core_counts = defaultdict(int)
    seq_pos = -1
    # print(orient)
    for i in range(len(seq)):
        if seq[i] == '-':
            continue
        seq_pos += 1

        if i == aln_core_start:
            core_dist = 0
        elif i < aln_core_start:
            core_dist = len(attd[i:aln_core_start].replace('-', ''))
        elif i > aln_core_start:
            core_dist = len(attd[aln_core_start:i].replace('-', ''))

        if orient == '3p':

            endseq_dist = len(nogaps_seq) - seq_pos
            start_genomic = core_dist - endseq_dist + pos + 1
            end_genomic = start_genomic + core_end - core_start
            cut_dist = start_genomic - pos - 1

        elif orient == '5p':
            # print("POS:", pos, "CORE_DIST:", core_dist, "SEQ_POS:", seq_pos)
            start_genomic = pos - core_dist + seq_pos
            end_genomic = start_genomic + core_end - core_start
            cut_dist = pos - end_genomic

        genomic_core_counts[(start_genomic, end_genomic, cut_dist)] += 1

    # print(genomic_core_counts)
    start_genomic, end_genomic, cut_dist = max(list(genomic_core_counts.items()), key=lambda x: x[1])[0]

    return start_genomic, end_genomic, cut_dist

def get_genomic_core(s, pos, orient, core_start, core_end, attd_fwd, attd_rev):

    matching_kmers = set()
    genomic_core_counts = defaultdict(int)
    primer_start_counts = defaultdict(int)

    pident_trimmed, pident_contig, align_length_trimmed, align_length_contig, trimmed_start, trimmed_end, contig_start, contig_end = None, None, None, None, None, None, None, None

    align1 = (None, None)
    if orient == '5p':
        attd = attd_rev
        core_start, core_end = len(attd) - core_end, len(attd) - core_start
        for p, kmer in enumerate(get_kmers(s, k=9)):

            attd_pos = attd.find(kmer)
            if attd_pos == -1:
                continue
            matching_kmers.add(kmer)
            start_coredist = attd_pos - core_start
            end_coredist = attd_pos - core_end
            start_genomic = pos - start_coredist + p
            end_genomic = pos - end_coredist + p
            cut_dist = pos - end_genomic
            genomic_core_counts[(start_genomic, end_genomic, cut_dist)] += 1

            primer_start = p + len(attd) - attd_pos
            primer_start_counts[primer_start] += 1

        if len(matching_kmers) > 0:

            primer_start = max(list(primer_start_counts.items()), key=lambda x: x[1])[0]
            if primer_start < len(s) and primer_start > 0:
                s_trimmed = s[:primer_start]
            else:
                s_trimmed = s

            align1 = \
            pairwise2.align.localms(s_trimmed, attd, match=2, mismatch=-1, open=-2, extend=-2, penalize_end_gaps=False)[
                0]
            pident_trimmed, align_length_trimmed, pident_contig, align_length_contig, trimmed_start, trimmed_end, contig_start, contig_end = get_percent_ident(
                align1[0], align1[1])
            trimmed_start, trimmed_end = len(attd) - trimmed_end, len(attd) - trimmed_start
            contig_start, contig_end = len(attd) - contig_end, len(attd) - contig_start
    else:
        attd = attd_fwd
        for p, kmer in enumerate(get_kmers(s, k=9)):
            attd_pos = attd.find(kmer)
            if attd_pos == -1:
                continue

            matching_kmers.add(kmer)
            start_coredist = core_start - attd_pos
            end_coredist = core_end - attd_pos
            endseq_dist = len(s) - p

            start_genomic = start_coredist - endseq_dist + pos + 1
            end_genomic = end_coredist - endseq_dist + pos + 1

            cut_dist = start_genomic - pos - 1
            genomic_core_counts[(start_genomic, end_genomic, cut_dist)] += 1

            primer_start = p - attd_pos
            primer_start_counts[primer_start] += 1

        if len(matching_kmers) > 0:
            primer_start = max(list(primer_start_counts.items()), key=lambda x: x[1])[0]
            if primer_start < len(s) and primer_start > 0:
                s_trimmed = s[primer_start:]
            else:
                s_trimmed = s

            align1 = \
            pairwise2.align.localms(s_trimmed, attd, match=2, mismatch=-1, open=-2, extend=-2, penalize_end_gaps=False)[
                0]
            pident_trimmed, align_length_trimmed, pident_contig, align_length_contig, trimmed_start, trimmed_end, contig_start, contig_end = get_percent_ident(
                align1[0], align1[1])
    try:
        start_genomic, end_genomic, cut_dist = max(list(genomic_core_counts.items()), key=lambda x: x[1])[0]
    except:
        return None

    align_start_genomic, align_end_genomic, align_cut_dist = get_align_core(align1[0], align1[1], pos, orient,
                                                                            core_start, core_end)

    if start_genomic != align_start_genomic:
        print(orient, start_genomic, end_genomic, cut_dist, align_start_genomic, align_end_genomic, align_cut_dist, sep
              ='\t')
    return start_genomic, end_genomic, cut_dist, align_start_genomic, align_end_genomic, align_cut_dist, len(
        matching_kmers), pident_trimmed, pident_contig, align_length_trimmed, align_length_contig, trimmed_start, trimmed_end, contig_start, contig_end, ';'.join(
        sorted(list(matching_kmers))), align1[0], align1[1]

def kmer_flanks(input, output, params):
    outbed = open(output.bed, 'w')
    outtsv = open(output.tsv, 'w')

    find = pd.read_csv(input.find, sep='\t')

    primer_start = None
    core_start, core_end = None, None
    attd_start, attd_end = None, None
    with open(input.donor_bed) as infile:
        for line in infile:
            line = line.strip().split()
            if line[3] == "Primer":
                primer_start = int(line[1])
            elif line[3] == 'Core':
                core_start, core_end = int(line[1]), int(line[2])
            elif line[3] == 'attD':
                attd_start, attd_end = int(line[1]), int(line[2])

    attd_fwd = None
    for rec in SeqIO.parse(input.donor_fasta, 'fasta'):
        attd_fwd = str(rec.seq)[primer_start:attd_end + 20]

    core_start, core_end = core_start - primer_start, core_end - primer_start
    # print(core_start, core_end)
    attd_rev = revcomp(attd_fwd)

    print('sample_id', 'name', 'softclip_pos', 'orient', 'softclip_count', 'chrom', 'kmer_core_start', 'kmer_core_end',
          'kmer_cut_dist', 'align_core_start', 'align_core_end', 'align_cut_dist', 'num_matching_kmers',
          'pident_trimmed', 'pident_contig', 'align_length_trimmed', 'align_length_contig', 'start_trimmed',
          'end_trimmed', 'start_contig', 'end_contig', 'matching_kmers', 'align_clipped', 'align_attd', sep
          ='\t', file=outtsv)

    for c, s, p, o, c5p, c3p in zip(find.contig, find.consensus_seq, find.pos, find.orient, find.softclip_count_5p,
                                    find.softclip_count_3p):
        name = c + ':' + str(p) + '(' + o + ')'
        read_count = c5p
        if c3p > c5p:
            read_count = c3p

        core = get_genomic_core(s, p, o, core_start, core_end, attd_fwd, attd_rev)
        if core is None:
            continue
        kmer_start_genomic, kmer_end_genomic, kmer_cut_dist, align_start_genomic, align_end_genomic, align_cut_dist, num_matching_kmers, pident_trimmed, pident_contig, align_length_trimmed, align_length_contig, trimmed_start, trimmed_end, contig_start, contig_end, matching_kmers, align_clipped, align_attd = core
        strand = '+'
        if o == '3p':
            strand = '-'

        print(c, kmer_start_genomic, kmer_end_genomic,
              'Name=' + name + ';method=kmer;Run=' + params.sample + ';kmer_cut_dist=' + str(
                  kmer_cut_dist) + ';' + 'read_count=' + str(read_count) + ';', kmer_cut_dist, strand, sep
              ='\t', file=outbed)
        print(c, align_start_genomic, align_end_genomic,
              'Name=' + name + ';method=align;Run=' + params.sample + ';align_cut_dist=' + str(
                  align_cut_dist) + ';' + 'read_count=' + str(read_count) + ';', align_cut_dist, strand, sep
              ='\t', file=outbed)
        print(params.sample, name, p, strand, read_count, c, kmer_start_genomic, kmer_end_genomic, kmer_cut_dist,
              align_start_genomic, align_end_genomic, align_cut_dist, num_matching_kmers, pident_trimmed, pident_contig,
              align_length_trimmed, align_length_contig, trimmed_start, trimmed_end, contig_start, contig_end,
              matching_kmers, align_clipped, align_attd, sep='\t', file=outtsv)
    outbed.close()
    outtsv.close()

    tsv = pd.read_csv(output.tsv, sep='\t', low_memory=False)
    tsv_tmp = tsv[['chrom', 'align_core_start', 'align_core_end', 'orient', 'softclip_count']]
    tsv_tmp['name'] = tsv_tmp['chrom'] + ':' + tsv_tmp['align_core_start'].astype(str) + '-' + tsv_tmp[
        'align_core_end'].astype(str) + '(' + tsv_tmp['orient'] + ')'
    tsv_tmp = tsv_tmp[['chrom', 'align_core_start', 'align_core_end', 'name', 'softclip_count', 'orient']]
    bed = BedTool.from_dataframe(tsv_tmp).slop(b=500, g=input.genome_file).sort()
    merged = bed.merge()
    final = bed.intersect(merged, wa=True, wb=True)
    final = pd.read_table(final.fn, names=['chrom', 'align_core_start', 'align_core_end', 'name', 'score', 'strand',
                                           'merged_chrom', 'merged_start', 'merged_end'])

    final['locus'] = final['merged_chrom'] + ':' + final['merged_start'].astype(str) + '-' + final['merged_end'].astype(
        str)

    locus_map = final[['name', 'locus']].drop_duplicates()

    tsv_tmp = tsv_tmp.merge(locus_map, on='name')
    outtsv = tsv.merge(tsv_tmp[['chrom', 'align_core_start', 'align_core_end', 'locus']].drop_duplicates(),
                       on=['chrom', 'align_core_start', 'align_core_end'])

    outtsv.to_csv(output.tsv, sep='\t', index=False)

def count_umis(input, output, params):
    tsv = pd.read_csv(input.tsv, sep='\t')
    loci = list(tsv['locus'].drop_duplicates())
    outbed = output.tsv + '.loci.bed'
    out = open(outbed, 'w')
    for l in loci:
        chrom = l.split(':')[0]
        start, end = l.split(':')[-1].split('-')
        print(chrom, start, end, l, sep='\t', file=out)
    out.close()

    outcounts = output.tsv + '.loci.umi_counts.tsv'
    cmd = "cat {bed} | bedtools sort -g {gen} | bedtools intersect -sorted -a stdin -b {bam} -g {gen} -wa -wb | cut -f4,8 | sed -r 's/\/.+//g' | sort | uniq | sed -r 's/\t.+UMI_/\t/g' | sort | uniq -c | sed -r 's/^ +//g' | sed -r 's/ /\t/g' > {out}".format(
        bed=outbed, gen=params.genome_file, bam=input.bam, out=outcounts)
    shell(cmd)

    try:
        counts = pd.read_csv(outcounts, sep='\t', header=None)
        counts.columns = ['count', 'locus', 'umi']
    except:
        counts = pd.DataFrame(columns=['count', 'locus', 'umi'])

    umi_counts = defaultdict(dict)
    for count, locus, umi in counts.values.tolist():
        if umi.count('N') > 1: continue
        if len(umi) != 12: continue
        umi_counts[locus][umi.encode('utf-8')] = count

    clusterer = UMIClusterer(cluster_method="directional")

    out = open(output.tsv, 'w')
    print("locus", "umi_counts", sep='\t', file=out)
    for locus in umi_counts:
        umis = umi_counts[locus]

        clustered_umis = clusterer(umis, threshold=1)

        print(locus, len(clustered_umis), sep='\t', file=out)
    out.close()

    os.remove(outcounts)
    os.remove(outbed)

def donor_count_umis(input, output, params):
    umi_counts = defaultdict(lambda: defaultdict(int))
    raw_counts = defaultdict(int)
    for read in pysam.AlignmentFile(input.bam, 'rb'):
        if read.is_reverse:
            continue
        if not read.is_paired:
            continue
        if read.tlen < 0:
            continue
        if read.tlen < 100:
            continue
        umi = read.query_name.split('UMI_')[-1]
        if len(umi) != 12: continue
        umi_counts[(read.reference_start, read.tlen)][umi.encode('utf-8')] += 1
        raw_counts[(read.reference_start, read.tlen)] += 1

    outfile = open(output.tsv, 'w')
    print("sample", "donor_start", "tlen", "UMI", "RAW", sep='\t', file=outfile)
    clusterer = UMIClusterer(cluster_method="directional")
    for locus in umi_counts:
        umis = umi_counts[locus]
        umis = dict(umis)
        clustered_umis = clusterer(umis, threshold=1)
        print(params.sample, *locus, len(clustered_umis), raw_counts[locus], sep='\t', file=outfile)
    outfile.close()

def count_raw(input, output, params):
    tsv = pd.read_csv(input.tsv, sep='\t')
    loci = list(tsv['locus'].drop_duplicates())
    outbed = output.tsv + '.loci.bed'
    out = open(outbed, 'w')
    for l in loci:
        chrom = l.split(':')[0]
        start, end = l.split(':')[-1].split('-')
        print(chrom, start, end, l, sep='\t', file=out)
    out.close()

    outcounts1 = output.tsv + '.loci.raw_counts.tsv'
    cmd = "cat {bed} | bedtools sort -g {gen} | bedtools intersect -sorted -a stdin -b {bam} -g {gen} -wa -wb | cut -f4,8 | sed -r 's/\/.+//g' | sort | uniq | cut -f1 | uniq -c | sort -g -k1 | sed -r 's/^ +//g' | sed -r 's/ /\t/g' > {out}".format(
        bed=outbed, gen=params.genome_file, bam=input.bam, out=outcounts1)
    shell(cmd)

    outcounts2 = output.tsv + '.loci.uniq_counts.tsv'
    cmd = "cat {bed} | bedtools sort -g {gen} | bedtools intersect -sorted -a stdin -b {bam} -g {gen} -wa -wb | cut -f4,8 | sed -r 's/\/.+//g' | sort | uniq | cut -f1 | uniq -c | sort -g -k1 | sed -r 's/^ +//g' | sed -r 's/ /\t/g' > {out}".format(
        bed=outbed, gen=params.genome_file, bam=input.rmdup_bam, out=outcounts2)
    shell(cmd)

    try:
        raw_counts = pd.read_csv(outcounts1, sep='\t', header=None)
        raw_counts.columns = ['count', 'locus']
    except:
        raw_counts = pd.DataFrame(columns=['count', 'locus'])

    try:
        uniq_counts = pd.read_csv(outcounts2, sep='\t', header=None)
        uniq_counts.columns = ['count', 'locus']
    except:
        uniq_counts = pd.DataFrame(columns=['count', 'locus'])

    raw_counts['count_type'] = 'RAW'
    uniq_counts['count_type'] = 'UNIQ'

    counts = pd.concat([raw_counts, uniq_counts])
    counts = counts[['count_type', 'locus', 'count']]

    counts.to_csv(output.tsv, sep='\t', index=None)

    os.remove(outcounts1)
    os.remove(outcounts2)

def get_integration_sites(input, output, params):
    cmd = 'cat {} | sed -r "s/Name=[^\t]+\t[^\t]+/Name\t0/g" | sort | uniq | bedtools sort | bedtools slop -b {} -g {} | bedtools getfasta -s -bed stdin -fi {} > {}'.format(
        input.bed, params.slop, input.genome_file, input.genome, output.sites)
    shell(cmd)

    tsv = pd.read_csv(input.kmer_flanks, sep='\t')
    counts = pd.read_csv(input.raw_counts, sep='\t')
    counts = counts[counts.count_type == 'UNIQ']
    counts = counts[['locus', 'count']]
    counts.columns = ['locus', 'uniq_counts']

    tsv = tsv.merge(counts, on='locus')
    tsv = tsv[(tsv.align_cut_dist >= -15) & (tsv.align_cut_dist <= 15)]
    tsv = tsv[(tsv.pident_trimmed > 0.8) | (tsv.pident_contig > 0.8)]
    tsv = tsv.sort_values(['uniq_counts', 'softclip_count'], ascending=[False, False])

    loci = list()
    for locus in tsv.locus:
        if locus not in loci:
            loci.append(locus)

    keep = set()
    keep_top = set()
    for locus in loci:
        tsv_filt = tsv[tsv.locus == locus]
        tsv_filt = tsv_filt.iloc[0, :]
        name = tsv_filt.chrom + ':' + str(tsv_filt.align_core_start - params.slop) + '-' + str(
            tsv_filt.align_core_end + params.slop) + '(' + tsv_filt.orient + ')'
        keep.add(name)
        if len(keep_top) < params.top:
            keep_top.add(name)

    out1 = open(output.curated, 'w')
    out2 = open(output.top, 'w')
    for rec in SeqIO.parse(output.sites, 'fasta'):
        if rec.id in keep:
            print('>' + str(rec.id), str(rec.seq), sep='\n', file=out1)
        if rec.id in keep_top:
            print('>' + str(rec.id), str(rec.seq), sep='\n', file=out2)
    out1.close()
    out2.close()

def get_site_annotations(input, output, params):
    bedfile = output.annot + '.tmp.bed'
    shell("cat {} | grep -v locus$ | cut -f6,10,11 | bedtools sort -g {} > {}".format(input.kmer_flanks,
                                                                                      params.genome_file, bedfile))

    gene_overlap = output.annot + '.tmp.gene_overlap.bed'
    shell("bedtools intersect -g {} -sorted -a {} -b {} -wa -wb > {}".format(params.genome_file, bedfile, params.genes,
                                                                             gene_overlap))

    exon_overlap = output.annot + '.tmp.exon_overlap.bed'
    shell("bedtools intersect -g {} -sorted -a {} -b {} -wa -wb > {}".format(params.genome_file, bedfile, params.exons,
                                                                             exon_overlap))

    overlaps_gene = defaultdict(set)
    with open(gene_overlap) as infile:
        for line in infile:
            line = line.strip().split('\t')
            sid = tuple(line[:3])
            gid = line[-1].split('gene_id=')[-1].split(';')[0]
            overlaps_gene[sid].add(gid)

    overlaps_exon = defaultdict(set)
    with open(exon_overlap) as infile:
        for line in infile:
            line = line.strip().split('\t')
            sid = tuple(line[:3])
            gid = line[-1].split('gene_id=')[-1].split(';')[0]
            overlaps_exon[sid].add(gid)

    annots = dict()

    for sid in overlaps_gene:
        if sid not in overlaps_exon:
            annots[sid] = ('INTRONIC', ';'.join(sorted(list(overlaps_gene[sid]))))
            continue
        exon_genes = set()
        for gid in overlaps_gene[sid]:
            if gid in overlaps_exon[sid]:
                exon_genes.add(gid)
        annots[sid] = ('EXONIC', ';'.join(sorted(list(exon_genes))))

    out = open(output.annot, 'w')
    print('chrom', 'start', 'end', 'region', 'gene', sep='\t', file=out)
    with open(bedfile) as infile:
        for line in infile:
            line = line.strip().split('\t')
            sid = tuple(line)
            if sid in annots:
                print(*sid, *annots[sid], sep='\t', file=out)
            else:
                if len(sid[0]) <= 5:
                    print(*sid, "INTERGENIC", "NA", sep='\t', file=out)
                else:
                    print(*sid, "UNANNOTATED", "NA", sep='\t', file=out)
    out.close()

    os.remove(bedfile)
    os.remove(exon_overlap)
    os.remove(gene_overlap)

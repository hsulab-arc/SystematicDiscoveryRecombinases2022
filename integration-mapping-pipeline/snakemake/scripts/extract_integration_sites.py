import subprocess
import os
import sys

OUTDIR = sys.argv[1]
INPATH = sys.argv[2]
GENOME = sys.argv[3]
OUTFILE = sys.argv[4]

sites_path = os.path.join(OUTDIR, "all_sites.bed")

outfile = open(sites_path, "w")
with open(INPATH) as infile:
    infile.readline()
    for line in infile:
        site = line.strip().split("\t")[5]
        chrom = site.split(":")[0]
        start, end = site.split(":")[1].split("(")[0].split("-")
        start, end = int(start) - 50, int(end) + 50
        strand = site.split("(")[-1].strip(")")
        print(chrom, start, end, site.split("(")[0], 0, strand, sep="\t", file=outfile)
outfile.close()

cmd = "sort {sites} | uniq | bedtools sort | bedtools getfasta -fi {genome} -bed stdin -s -nameOnly > {outfile}".format(
    sites=sites_path, genome=GENOME, outfile=OUTFILE
)
subprocess.run(cmd, shell=True)
os.remove(sites_path)

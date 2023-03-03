require(tidyverse)
require(bedr)
require(HelloRanges)
require(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)
WORKDIR = args[1]
#WORKDIR <- "/home/alisonf/uditas-pipeline/WORKDIR_pool2/"
OUTDIR = args[2]
#OUTDIR <- "/home/alisonf/uditas-pipeline/WORKDIR_pool2/11.results"

# Load the metadata
metadata <- read_tsv(paste0(WORKDIR, '/metadata.tsv')) %>%
  dplyr::select(sample_id=sample, lsr, biorep, techrep) %>%
  mutate(biorep=paste0('biorep', biorep), techrep=paste0('techrep', techrep))

# Get raw counts
all_counts <- NULL
for (f in Sys.glob(paste0(WORKDIR, "/09.raw_counts/*tsv"))){
  samp <- gsub('\\..+', '', basename(f))
  df <- read_tsv(f) %>% mutate(sample_id=samp) %>%
    dplyr::select(sample_id, count_type, locus, count)
  if (is.null(all_counts)){
    all_counts <- df
  }else{
    all_counts <- rbind(all_counts, df)
  }
}
for (f in Sys.glob(paste0(WORKDIR, "/09.umi_counts/*.umi_counts.tsv"))){
  samp <- gsub('\\..+', '', basename(f))
  df <- read_tsv(f) %>% mutate(sample_id=samp, count_type='UMI') %>%
    dplyr::select(sample_id, count_type, locus, count=umi_counts)
  if (is.null(all_counts)){
    all_counts <- df
  }else{
    all_counts <- rbind(all_counts, df)
  }
}

flanks <- NULL
for (f in Sys.glob(paste0(WORKDIR, "/08.kmer_flanks/*kmer_flanks.tsv"))){
  samp <- gsub('\\..+', '', basename(f))
  df <- read_tsv(f) %>% mutate(sample_id=samp) %>%
    dplyr::select(sample_id, name, locus, chrom, start=align_core_start,
                  end=align_core_end, orient, find_count=softclip_count,
                  align_cut_dist, pident_trimmed, pident_contig,
                  align_length_trimmed, align_length_contig, start_trimmed,
                  end_trimmed, start_contig, end_contig)
  if (is.null(flanks)){
    flanks <- df
  }else{
    flanks <- rbind(flanks, df)
  }
}

top_sites <- flanks %>% group_by(sample_id, locus) %>%
  filter(find_count==max(find_count)) %>%
  dplyr::count(sample_id, locus, chrom, start, end, orient) %>% filter(n==max(n)) %>%
  filter(row_number() == 1) %>% group_by() %>%
  mutate(site=paste(chrom, paste(start, end, sep='-'), orient, sep=':')) %>%
  dplyr::select(sample_id, locus, site)

keep_loci <- flanks %>%
  filter(align_cut_dist >= -15, align_cut_dist <= 15) %>%
  filter(pident_trimmed > 0.8 | pident_contig > 0.8) %>%
  dplyr::select(sample_id, locus) %>% unique()

results <- all_counts %>%
  spread(count_type, count, fill=0) %>%
  arrange(locus)

unfiltered_counts <- results %>% inner_join(top_sites) %>%
  group_by(sample_id) %>% arrange(desc(UNIQ)) %>%
  mutate(rank_uniq_reads=row_number()) %>%
  mutate(perc_uniq_reads=UNIQ/sum(UNIQ)*100) %>%
  arrange(desc(RAW)) %>%
  mutate(rank_raw_reads=row_number()) %>%
  mutate(perc_raw_reads=RAW/sum(RAW)*100) %>%
  arrange(desc(UMI)) %>%
  mutate(rank_umi_reads=row_number()) %>%
  mutate(perc_umi_reads=UMI/sum(UMI)*100) %>%
  group_by()

filtered_counts <- results %>% inner_join(keep_loci) %>%
  inner_join(top_sites) %>% inner_join(keep_loci) %>%
  group_by(sample_id) %>%
  arrange(desc(UNIQ)) %>%
  mutate(rank_uniq_reads=row_number()) %>%
  mutate(perc_uniq_reads=UNIQ/sum(UNIQ)*100) %>%
  arrange(desc(RAW)) %>%
  mutate(rank_raw_reads=row_number()) %>%
  mutate(perc_raw_reads=RAW/sum(RAW)*100) %>%
  arrange(desc(UMI)) %>%
  mutate(rank_umi_reads=row_number()) %>%
  mutate(perc_umi_reads=UMI/sum(UMI)*100) %>%
  group_by()


unfiltered_counts <- unfiltered_counts %>% mutate(site=gsub(':(.)$', '(\\1)', site))
filtered_counts <- filtered_counts %>% mutate(site=gsub(':(.)$', '(\\1)', site))
write_tsv(unfiltered_counts, paste0(OUTDIR, "/raw_counts.unfiltered.tsv"))
write_tsv(filtered_counts, paste0(OUTDIR, "/raw_counts.filtered.tsv"))

# Get site_annotations
annots <- NULL
for (f in Sys.glob(paste0(WORKDIR, "/10.target_sites/*annotations.tsv"))){
  samp <- gsub('\\..+', '', basename(f))
  df <- read_tsv(f) #%>% mutate(sample_id=samp) %>%
    #dplyr::select(sample_id, chrom, start, end, region, gene)
  if (ncol(df) == 0){
  	df <- data.frame(sample_id=character(), chrom=character(), start=integer(), end=integer(), region=character(), gene=character())
  }else{
        df <- df %>% mutate(sample_id=samp) %>%
		dplyr::select(sample_id, chrom, start, end, region, gene)
  }
  if (is.null(annots)){
    annots <- df
  }else{
    annots <- rbind(annots, df)
  }
}

annots <- rbind(
  annots %>% mutate(site=paste0(chrom, ':', start, '-', end, '(+)')),
  annots %>% mutate(site=paste0(chrom, ':', start, '-', end, '(-)'))
) %>% dplyr::select(sample_id, chrom, start, end, site, region, gene) %>%
  arrange(site) %>%
  inner_join(unfiltered_counts %>% dplyr::select(site) %>% unique())
annots %>% write_tsv(paste0(OUTDIR, "/site_annotations.tsv"))


# Get donor-specific UMI count
donor_umi_counts <- NULL
for (f in Sys.glob(paste0(WORKDIR, "/09.umi_counts/*donor*tsv"))){
  samp <- gsub('\\..+', '', basename(f))
  df <- read_tsv(f) %>% dplyr::rename(sample_id=sample) %>%
    dplyr::select(sample_id, donor_start, tlen, donor_UMI=UMI, donor_RAW=RAW)
  if (is.null(donor_umi_counts)){
    donor_umi_counts <- df
  }else{
    donor_umi_counts <- rbind(donor_umi_counts, df)
  }
}
donor_umi_counts <- donor_umi_counts %>% group_by(sample_id) %>%
  summarize(donor_UMI=sum(donor_UMI), donor_RAW=sum(donor_RAW), donor_UNIQ=n())

donor_umi_counts %>% write_tsv(paste0(OUTDIR, "/donor_umi_counts.tsv"))

# Now merge the counts
unique_loci <- filtered_counts %>% dplyr::select(locus) %>% unique() %>%
  group_by()
tmp <- unique_loci %>% mutate(locus2=locus) %>%
  separate(locus, into=c('chrom', 'start', 'end'), sep='[:-]') %>%
  mutate(start=as.numeric(start), end=as.numeric(end)) %>%
  unique()
ir <- IRanges(start=tmp$start, end=tmp$end, names=tmp$locus2)
gr <- GRanges(tmp$chrom, ir)
ovrlps <- findOverlapPairs(GenomicRanges::reduce(gr), gr)
ovrlps <- data.frame(
  merged_chrom=as.character(ovrlps@first@seqnames),
  merged_start=ovrlps@first@ranges@start,
  merged_end=ovrlps@first@ranges@start + ovrlps@first@ranges@width - 1,
  chrom=as.character(ovrlps@second@seqnames),
  start=ovrlps@second@ranges@start,
  end=ovrlps@second@ranges@start + ovrlps@second@ranges@width - 1
) %>% group_by()
merged_loci <- tmp %>% inner_join(ovrlps) %>%
  mutate(merged_locus=paste(merged_chrom, paste(merged_start, merged_end, sep='-'), sep=':')) %>%
  dplyr::select(merged_locus, locus=locus2) %>%
  dplyr::select(merged_locus, locus)

merged_sample <- filtered_counts %>% inner_join(merged_loci) %>% inner_join(metadata) %>%
  dplyr::select(sample_id, lsr, biorep, techrep, locus, merged_locus, site,
                RAW, UNIQ, UMI) %>%
  arrange(desc(UMI)) %>% group_by(sample_id) %>%
  mutate(perc_umi=UMI/sum(UMI)*100) %>%
  mutate(umi_rank=row_number()) %>% group_by()
write_tsv(merged_sample, paste0(OUTDIR, "/merged_counts.sample.tsv"))

merged_biorep <- merged_sample %>% group_by(lsr, biorep, merged_locus) %>%
  summarize(site=site[UMI==max(UMI)][1], RAW=sum(RAW), UNIQ=sum(UNIQ), UMI=sum(UMI)) %>%
  group_by() %>% arrange(desc(UMI)) %>% group_by(lsr, biorep) %>%
  mutate(perc_umi=UMI/sum(UMI)*100) %>%
  mutate(umi_rank=row_number()) %>% group_by()
write_tsv(merged_biorep, paste0(OUTDIR, "/merged_counts.biorep.tsv"))

merged_lsr <- merged_biorep %>% group_by(lsr, merged_locus) %>%
  summarize(site=site[UMI==max(UMI)][1], RAW=sum(RAW), UNIQ=sum(UNIQ), UMI=sum(UMI)) %>%
  group_by() %>% arrange(desc(UMI)) %>% group_by(lsr) %>%
  mutate(perc_umi=UMI/sum(UMI)*100) %>%
  mutate(umi_rank=row_number()) %>% group_by()
write_tsv(merged_lsr, paste0(OUTDIR, "/merged_counts.lsr.tsv"))

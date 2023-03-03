require(tidyverse)
require(Biostrings)
require(ggseqlogo)

args = commandArgs(trailingOnly=TRUE)
OUTDIR = args[1]
#OUTDIR <- "/home/mdurrant/code/uditas-pipeline/test/WORKDIR/11.results"
MERGED_SAMPLE = args[2]
#MERGED_SAMPLE <- "/home/mdurrant/code/uditas-pipeline/test/WORKDIR/11.results/merged_counts.sample.tsv"
MERGED_BIOREP = args[3]
#MERGED_BIOREP <- "/home/mdurrant/code/uditas-pipeline/test/WORKDIR/11.results/merged_counts.biorep.tsv"
MERGED_LSR = args[4]
#MERGED_LSR <- "/home/mdurrant/code/uditas-pipeline/test/WORKDIR/11.results/merged_counts.lsr.tsv"
SITE_SEQS = args[5]
#SITE_SEQS <- "/home/mdurrant/code/uditas-pipeline/test/WORKDIR/11.results/all_integration_sites.fna"

write(paste(OUTDIR, MERGED_SAMPLE, MERGED_BIOREP, MERGED_LSR, SITE_SEQS, sep='\n'), stdout())
siteseqs <- data.frame(site=names(readAAStringSet(SITE_SEQS)),
           seq=as.character(readAAStringSet(SITE_SEQS))) %>% 
  group_by() %>% mutate(seq=toupper(seq))

merged_sample <- read_tsv(MERGED_SAMPLE) %>% 
  inner_join(siteseqs)

outdir_merged_sample <- paste0(OUTDIR, '/seqlogos/sample') 
dir.create(outdir_merged_sample, showWarnings=F, recursive=T)
for (samp in unique(merged_sample$sample_id)){
  filt <- merged_sample %>% filter(sample_id==samp)
  
  top_20 <- filt %>% filter(umi_rank <= 20) %>% filter(!grepl('[^ACGT]', seq)) %>% 
    .$seq
  gg_20 <- ggseqlogo(top_20) + theme_classic() +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 2)) +
    ggtitle(paste0(samp, '; n=', length(top_20)))
  
  top_100 <- filt %>% filter(umi_rank <= 100) %>% filter(!grepl('[^ACGT]', seq)) %>% 
    .$seq
  gg_100 <- ggseqlogo(top_100) + theme_classic() +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 2)) +
    ggtitle(paste0(samp, '; n=', length(top_100)))
  
  top_200 <- filt %>% filter(umi_rank <= 200) %>% filter(!grepl('[^ACGT]', seq)) %>% 
    .$seq
  gg_200 <- ggseqlogo(top_200) + theme_classic() +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 2)) +
    ggtitle(paste0(samp, '; n=', length(top_200)))
  
  top_all <- filt %>% filter(!grepl('[^ACGT]', seq)) %>% .$seq
  gg_all <- ggseqlogo(top_all) + theme_classic() +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 2)) +
    ggtitle(paste0(samp, '; n=', length(top_all)))
  
  ggsave(paste0(outdir_merged_sample, '/', samp, '.top20.pdf'), plot=gg_20, units='mm', width=180, height=45)
  ggsave(paste0(outdir_merged_sample, '/', samp, '.top100.pdf'), plot=gg_100, units='mm', width=180, height=45)
  ggsave(paste0(outdir_merged_sample, '/', samp, '.top200.pdf'), plot=gg_200, units='mm', width=180, height=45)
  ggsave(paste0(outdir_merged_sample, '/', samp, '.all.pdf'), plot=gg_all, units='mm', width=180, height=45)
}

merged_biorep <- read_tsv(MERGED_BIOREP) %>% 
  inner_join(siteseqs) %>% mutate(biorep=paste(lsr, biorep, sep='-'))

outdir_merged_biorep <- paste0(OUTDIR, '/seqlogos/biorep')
dir.create(outdir_merged_biorep, showWarnings=F, recursive=T)
for (brep in unique(merged_biorep$biorep)){
  filt <- merged_biorep %>% filter(biorep==brep)
  
  top_20 <- filt %>% filter(umi_rank <= 20) %>% filter(!grepl('[^ACGT]', seq)) %>% 
    .$seq
  gg_20 <- ggseqlogo(top_20) + theme_classic() +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 2)) +
    ggtitle(paste0(brep, '; n=', length(top_20)))
  
  top_100 <- filt %>% filter(umi_rank <= 100) %>% filter(!grepl('[^ACGT]', seq)) %>% 
    .$seq
  gg_100 <- ggseqlogo(top_100) + theme_classic() +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 2)) +
    ggtitle(paste0(brep, '; n=', length(top_100)))
  
  top_200 <- filt %>% filter(umi_rank <= 200) %>% filter(!grepl('[^ACGT]', seq)) %>% 
    .$seq
  gg_200 <- ggseqlogo(top_200) + theme_classic() +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 2)) +
    ggtitle(paste0(brep, '; n=', length(top_200)))
  
  top_all <- filt %>% filter(!grepl('[^ACGT]', seq)) %>% .$seq
  gg_all <- ggseqlogo(top_all) + theme_classic() +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 2)) +
    ggtitle(paste0(brep, '; n=', length(top_all)))
  
  ggsave(paste0(outdir_merged_biorep, '/', brep, '.top20.pdf'), plot=gg_20, units='mm', width=180, height=45)
  ggsave(paste0(outdir_merged_biorep, '/', brep, '.top100.pdf'), plot=gg_100, units='mm', width=180, height=45)
  ggsave(paste0(outdir_merged_biorep, '/', brep, '.top200.pdf'), plot=gg_200, units='mm', width=180, height=45)
  ggsave(paste0(outdir_merged_biorep, '/', brep, '.all.pdf'), plot=gg_all, units='mm', width=180, height=45)
}

merged_lsr <- read_tsv(MERGED_LSR) %>% 
  inner_join(siteseqs)

outdir_merged_lsr <- paste0(OUTDIR, '/seqlogos/lsr')
dir.create(outdir_merged_lsr, showWarnings=F, recursive=T)
for (lsr_id in unique(merged_lsr$lsr)){
  filt <- merged_lsr %>% filter(lsr==lsr_id)
  
  top_20 <- filt %>% filter(umi_rank <= 20) %>% filter(!grepl('[^ACGT]', seq)) %>% 
    .$seq
  gg_20 <- ggseqlogo(top_20) + theme_classic() +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 2)) +
    ggtitle(paste0(lsr_id, '; n=', length(top_20)))
  
  top_100 <- filt %>% filter(umi_rank <= 100) %>% filter(!grepl('[^ACGT]', seq)) %>% 
    .$seq
  gg_100 <- ggseqlogo(top_100) + theme_classic() +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 2)) +
    ggtitle(paste0(lsr_id, '; n=', length(top_100)))
  
  top_200 <- filt %>% filter(umi_rank <= 200) %>% filter(!grepl('[^ACGT]', seq)) %>% 
    .$seq
  gg_200 <- ggseqlogo(top_200) + theme_classic() +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 2)) +
    ggtitle(paste0(lsr_id, '; n=', length(top_200)))
  
  top_all <- filt %>% filter(!grepl('[^ACGT]', seq)) %>% .$seq
  gg_all <- ggseqlogo(top_all) + theme_classic() +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 2)) +
    ggtitle(paste0(lsr_id, '; n=', length(top_all)))
  
  ggsave(paste0(outdir_merged_lsr, '/', lsr_id, '.top20.pdf'), plot=gg_20, units='mm', width=180, height=45)
  ggsave(paste0(outdir_merged_lsr, '/', lsr_id, '.top100.pdf'), plot=gg_100, units='mm', width=180, height=45)
  ggsave(paste0(outdir_merged_lsr, '/', lsr_id, '.top200.pdf'), plot=gg_200, units='mm', width=180, height=45)
  ggsave(paste0(outdir_merged_lsr, '/', lsr_id, '.all.pdf'), plot=gg_all, units='mm', width=180, height=45)
}

write("All seqlogos completed successfully and can be found in 11.results/seqlogos.",
      paste0(OUTDIR, '/seqlogos.txt'))

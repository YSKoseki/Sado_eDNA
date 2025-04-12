# 02_DADA2_Lib2019.R
# An R script to perform DADA2 error-correction algorithm
# Reference https://benjjneb.github.io/dada2/index.html
# R 4.4.2

# Packages required
library(dada2); packageVersion("dada2") # 1.34.0
library(svglite); packageVersion("svglite") # 2.1.3
library(scales); packageVersion("scales") # 1.3.0
library(tidyverse); packageVersion("tidyverse") # 2.0.0

# Directory/File Paths for:
## Primer-trimmed fastq files
path_cutadapt1 <- "01_Cutadapt_Lib2019/FASTQ_trim"
## Primer trimming stats
path_cutadapt2 <- "01_Cutadapt_Lib2019/summary.tsv"
## Output files
path_out <- "02_DADA2_Lib2019"

# Library ID used as part of the sequence identifiers in the output FASTA file
lib_id <- "Lib2019"

# FASTQ file identifiers for forward and reverse reads
fastq_str_fwds <- "R1_001_trim.fastq.gz"
fastq_str_revs <- "R2_001_trim.fastq.gz"

# FASTQ file identifiers before and after quality filtering and trimming
fastq_str_pre <- "trim"
fastq_str_post <- "qlt"

# Plot the read quality profiles before quality filtering and trimming
## Output directory
path_qprfl <- paste0(path_out, "/01_qlt_prfl")
dir.create(path_qprfl, recursive = TRUE)
## Name lists of input FASTQ files (zero-read samples not including)
fwds <- path_cutadapt1 %>% 
  list.files(pattern = fastq_str_fwds, full.names = TRUE) %>%
  sort()
revs <- path_cutadapt1 %>% 
  list.files(pattern = fastq_str_revs, full.names = TRUE) %>%
  sort()
## Retrieve sample names - may change depending on the directory and file names
sample_names <- fwds %>% 
  strsplit("/") %>% map_vec(`[`, 3) %>% strsplit("_") %>% map_vec(`[`, 1)
## Read-containing samples
primtrim_tab <- read_tsv(path_cutadapt2)
zero_read_smpl <- primtrim_tab %>% 
  filter(totalRP == 0) %>% 
  select(sample) %>% 
  unlist()
is_zero <- sample_names %in% zero_read_smpl
nzr_sample_names <- sample_names[!is_zero]
nzr_fwds <- fwds[!is_zero]
nzr_revs <- revs[!is_zero]
## Plotting
map2(
  .x = nzr_sample_names, 
  .y = seq_along(nzr_sample_names), 
  .f = ~ {
    plotQualityProfile(c(nzr_fwds[.y], nzr_revs[.y]))
    ggsave(
      paste0(.x, ".svg"), 
      path = path_qprfl, 
      device = "svg", width = 20, height = 20, units = "cm"
    )
  }
)

# Quality filtering and trimming
## Output directory
path_cutadapt1 <- paste0(path_out, "/02_qlt_FASTQ")
dir.create(path_cutadapt1, recursive = TRUE)
## Name lists of output FASTQ files
nzr_fwds_flt <- nzr_fwds %>% 
  strsplit("/") %>% 
  map(`[`, 3) %>% 
  paste0(path_cutadapt1, "/", .) %>% 
  gsub(fastq_str_pre, fastq_str_post, .)
nzr_revs_flt <- nzr_revs %>% 
  strsplit("/") %>% 
  map(`[`, 3) %>% 
  paste0(path_cutadapt1, "/", .) %>% 
  gsub(fastq_str_pre, fastq_str_post, .)
## Do the quality filtering and trimming
## Set `minLen` to 90 bp: MiFish amplicons average ~172 bp and a >=12 bp overlap
##  needed for merging.
## Other settings remain at their defaults.
qfltrd <- filterAndTrim(
  nzr_fwds, nzr_fwds_flt, nzr_revs, nzr_revs_flt, 
  minLen = 90, truncQ = 2, maxN = 0, maxEE = c(2, 2), multithread = TRUE
)
## Simplify row names
rownames(qfltrd) <- nzr_sample_names

# Plot read quality profiles after quality filtering and trimming
## Output directory
path_qprfl2 <- paste0(path_out, "/03_qlt_prfl_flt/")
dir.create(path_qprfl2, recursive = TRUE)
## Identifier of zero-read samples, if any after quality filtering
zero_read_smpl2 <- which(qfltrd[, "reads.out"] == 0) %>% 
  names()
## Name lists of input FASTQ files (zero-read samples not including)
iszero2 <- nzr_sample_names %in% zero_read_smpl2
nzr_sample_names2 <- nzr_sample_names[!iszero2]
nzr_fwds_flt2 <- nzr_fwds_flt[!iszero2]
nzr_revs_flt2 <- nzr_revs_flt[!iszero2]
## Plotting
map2(
  .x = nzr_sample_names2, 
  .y = seq_along(nzr_sample_names2), 
  .f = ~ {
    plotQualityProfile(c(nzr_fwds_flt2[.y], nzr_revs_flt2[.y]))
    ggsave(
      paste0(.x, "_qlt.svg"), 
      path = path_qprfl2, 
      device = "svg", width = 20, height = 20, units = "cm"
    )
  }
)

# Generate base calling error models
## Learn error rates from FASTQ files
set.seed(101)  # for reproducibility
fwds_error <- nzr_fwds_flt2 %>% 
  learnErrors(
    multithread = TRUE, randomize = TRUE, verbose = TRUE, MAX_CONSIST = 10, 
    nbases = 1e+08 # All default values
  ) 
set.seed(101)
revs_error <- nzr_revs_flt2 %>% 
  learnErrors(
    multithread = TRUE, randomize = TRUE, verbose = TRUE, MAX_CONSIST = 10, 
    nbases = 1e+08 # All default values
  )
## Output directory
path_error_plot <- paste0(path_out, "/04_error_plot")
dir.create(path_error_plot, recursive = TRUE)
## Plot base calling error profiles
plotErrors(fwds_error, nominalQ = TRUE)
ggsave(
  "fwds_error.svg", path = path_error_plot, device = "svg", 
  width = 20, height = 20, units = "cm"
)
plotErrors(revs_error, nominalQ = TRUE)
ggsave(
  "revs_error.svg", path = path_error_plot, device = "svg", 
  width = 20, height = 20, units = "cm"
)

# Denoising
## Dereplicate sequencing reads
fwds_derep <- derepFastq(nzr_fwds_flt2, verbose = TRUE)
revs_derep <- derepFastq(nzr_revs_flt2, verbose = TRUE)
## Simplify row names
names(fwds_derep) <- names(revs_derep) <- nzr_sample_names2
## Denoise the dereplicated reads
set.seed(101)
fwds_dada <- dada(fwds_derep, err = fwds_error, multithread = TRUE)
set.seed(101)
revs_dada <- dada(revs_derep, err = revs_error, multithread = TRUE)
## Merge forward and reverse reads
mergers <- mergePairs(
  fwds_dada, fwds_derep, revs_dada, revs_derep, verbose = TRUE
)
## Sequence table
seq_tab <- makeSequenceTable(mergers)
## Sequence count by length
seq_tab %>% getSequences() %>% nchar() %>% table()

# Additional filtering
## Length filtering: Discard sequences outside the range of MiFish amplicon 
## lengths (160-190 bp)
seq_len <- seq_tab %>% colnames() %>% nchar()
seq_1en2 <- seq_len %in% 160:190
seq_tab2 <- seq_tab %>% `[`(, seq_1en2)
## Sequence count by length
seq_tab2 %>% getSequences() %>% nchar() %>% table()
## Remove chimeras
asv_tab <- seq_tab2 %>% 
  removeBimeraDenovo(multithread = TRUE, verbose = TRUE, method = "consensus")
## Sequence count by length
asv_tab %>% getSequences() %>% nchar() %>% table()

# Create a summary table of the Cutadapt-DADA2 pipeline
## Helper function to get the number of unique sequences
getN <- function(x) sum(getUniques(x))
## Summary table
summary_sub <- primtrim_tab %>% 
  mutate(sample = sample, initial = as.integer(totalRP), .keep = "none")
summary_sub2 <- tibble(
  sample = nzr_sample_names, 
  primtrim = qfltrd[, "reads.in"] %>% as.integer(),
  qltfilter = qfltrd[, "reads.out"] %>% as.integer()
)
summary_sub3 <- tibble(
  sample = nzr_sample_names2,
  dadafwds = map_vec(fwds_dada, getN), 
  dadarevs = map_vec(revs_dada, getN), 
  dadamerge = map_vec(mergers, getN), 
  lenfilter = rowSums(seq_tab2) %>% as.integer(), 
  chimremove = rowSums(asv_tab) %>% as.integer(), 
)
summary_tab <- summary_sub %>% 
  left_join(summary_sub2) %>% 
  left_join(summary_sub3) %>% 
  mutate(across(everything(), ~replace_na(., 0))) %>%
  mutate(pct_remain = percent(chimremove / initial, accuracy = 0.1))
summary_tab %>% 
  print(n = summary_sub %>% nrow())
write_tsv(summary_tab, file = paste0(path_out, "/DADA2_sum.tsv"))

# Save data and results
## Save ASV sequences as a FASTA file
asv_seqs <- colnames(asv_tab)
n_asvs <- length(asv_seqs)
asv_seqs_fa <- paste0(">", lib_id, "_%04d") %>% 
  sprintf(1:n_asvs) %>% 
  rbind(asv_seqs) %>% 
  c() %>% 
  as.matrix(ncol = 1)
asv_seqs_fa %>% 
  as.data.frame() %>% 
  write_tsv(file = paste0(path_out, "/ASV.fa"), col_names = FALSE)
## R objects
path_obj <- paste0(path_out, "/05_saved_obj")
dir.create(path_obj, recursive = TRUE)
saveRDS(asv_tab, paste0(path_obj, "/ASV_tab.obj"))
saveRDS(summary_tab, paste0(path_obj, "/summary_tab.obj"))
saveRDS(asv_seqs_fa, paste0(path_obj, "/ASV_fa.obj"))
## Workspace
paste0(path_out, "/saved.Rdata") %>% save.image()
## Session info
sessionInfo() %>% 
  capture.output() %>% 
  writeLines(paste0(path_out, "/session.info"))

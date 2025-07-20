# 05_Mergeruns.R
# R script to merge data from different sequencing runs
# R 4.2.2

# Packages required
library(dada2); packageVersion("dada2") # 1.34.0
library(phyloseq); packageVersion("phyloseq") # 1.50.0
library(Biostrings); packageVersion("Biostrings") # 2.74.1
library(tidyverse); packageVersion("tidyverse") # 2.0.0

# Paths
# Input data
path_in <- list(
  run2018i = "./01_run2018i/03_gmmDenoise_run2018i/ASV_tab_gmm.rds",
  run2018ii = "./02_run2018ii/03_gmmDenoise_run2018ii/ASV_tab_gmm.rds",
  run2019 = "./03_run2019/03_gmmDenoise_run2019/ASV_tab_gmm.rds",
  run2022 = "./04_run2022/03_gmmDenoise_run2022/ASV_tab_gmm.rds"
)
## Output directory
path_out <- "05_Mergeruns"

# Load R objects
ASV_tab_lis <- path_in %>% map(readRDS)

# Arrange ASV tables by year
## Custom function to label sample names with sampling years
smpl_yr <- function(tbl, yr) {
  newtbl <- tbl
  nam <- rownames(tbl)
  newnam <- paste0(yr, "-", nam)
  rownames(newtbl) <- newnam
  newtbl
}
## Organize by-year ASV tables into a nested tibble
ASV_tab_nest <- ASV_tab_lis[c("run2018i", "run2018ii")] %>% 
  dada2::mergeSequenceTables(tables = .) %>%
  list(
    run2018 = .,
    run2019 = ASV_tab_lis[["run2019"]],
    run2022 = ASV_tab_lis[["run2022"]]
  ) %>% 
  tibble(
    yr = names(.) %>% str_replace("run", "") %>% as.integer(),
    asvtbl = .
  ) %>% 
  mutate(
    asvtbl = map2(asvtbl, yr, ~ smpl_yr(.x, .y))
  )

# Merge ASV tables and convert the merged table into phyloseq objects
## Custom function to create a DNAStringSet object
add_seq <- function(ps, prefix = "ASV", sep = "", digit = 3) {
  nam <- ps %>% phyloseq::taxa_names()
  n <- length(nam)
  seq <- nam %>% Biostrings::DNAStringSet()
  names(seq) <- nam
  new_ps <- phyloseq::merge_phyloseq(ps, seq)
  phyloseq::taxa_names(new_ps) <- paste0(prefix, sep, "%0", digit, "d") %>% 
    sprintf(seq(n))
  new_ps
}
## Merge ASV tables
glob_ps <- ASV_tab_nest$asvtbl %>% 
  dada2::mergeSequenceTables(
    tables = ., repeats = "error", orderBy = "abundance"
  ) %>%
  phyloseq::otu_table(taxa_are_rows = FALSE) %>% 
  phyloseq() %>% 
  add_seq()

# Create FASTA file
make_fasta <- function(ps, filename = "ASV.fa") {
  seq <- ps %>% refseq() %>% as.character()
  nam <- names(seq)
  fasta <- paste0(">", nam) %>% 
    rbind(seq) %>% 
    c() %>% 
    as.matrix(ncol = 1)
  fasta %>% 
    as.data.frame() %>% 
    write_tsv(file = paste0(path_out, "/ASV.fa"), col_names = FALSE)
}
dir.create(path_out, recursive = TRUE)
make_fasta(glob_ps, filename = paste0(path_out, "/ASV.fa"))

# Outputs
## R objects
saveRDS(glob_ps, paste0(path_out, "/glob_ps.rds"))
## Workspace
paste0(path_out, "/saved.Rdata") %>% save.image()
## Session info
sessionInfo() %>% 
  capture.output() %>% 
  writeLines(paste0(path_out, "/session.info"))

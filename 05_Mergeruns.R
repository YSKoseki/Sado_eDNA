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
  run2018i = "./01_Run2018i/03_gmmDenoise_Run2018i/ASV_tab_gmm.rds",
  run2018ii = "./02_Run2018ii/03_gmmDenoise_Run2018ii/ASV_tab_gmm.rds",
  run2019 = "./03_Run2019/03_gmmDenoise_Run2019/ASV_tab_gmm.rds",
  run2022 = "./04_Run2022/03_gmmDenoise_Run2022/ASV_tab_gmm.rds",
  idmap = "./05_Mergeruns/ASV_id_map.tsv"
)
## Output directory
path_out <- "05_Mergeruns"

# Load R objects
ASV_tab_lis <- path_in[c("run2018i", "run2018ii", "run2019", "run2022")] %>%
  map(readRDS)

# Load the ASV ID map (sequence-to-ID table) from previous runs of this script,
#  so that existing sequences keep their IDs when new data are added
id_map <- if (file.exists(path_in$idmap)) {
  read_tsv(path_in$idmap, col_types = "cc")
} else {
  tibble(id = character(), seq = character())
}

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
## Custom function to create a DNAStringSet object and name ASVs; sequences
##  found in the ID map keep their IDs, and new sequences are numbered after
##  the largest existing ID in order of abundance
add_seq <- function(ps, id_map, prefix = "ASV", sep = "", digit = 3) {
  nam <- ps %>% phyloseq::taxa_names()
  seq <- nam %>% Biostrings::DNAStringSet()
  names(seq) <- nam
  new_ps <- phyloseq::merge_phyloseq(ps, seq)
  id <- id_map$id[match(nam, id_map$seq)]
  is_new <- is.na(id)
  max_id <- id_map$id %>%
    str_remove(paste0("^", prefix, sep)) %>%
    as.integer() %>%
    max(0)
  id[is_new] <- paste0(prefix, sep, "%0", digit, "d") %>%
    sprintf(max_id + seq_len(sum(is_new)))
  phyloseq::taxa_names(new_ps) <- id
  new_ps
}
## Merge ASV tables
glob_ps <- ASV_tab_nest$asvtbl %>%
  dada2::mergeSequenceTables(
    tables = ., repeats = "error", orderBy = "abundance"
  ) %>%
  phyloseq::otu_table(taxa_are_rows = FALSE) %>%
  phyloseq() %>%
  add_seq(id_map = id_map)
## Update the ID map, appending new sequences and keeping all previous entries
##  (including sequences absent from the current data) so IDs are never reused
id_map <- glob_ps %>%
  refseq() %>%
  as.character() %>%
  tibble(id = names(.), seq = .) %>%
  bind_rows(id_map, .) %>%
  distinct(seq, .keep_all = TRUE)

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
## ASV ID map
write_tsv(id_map, paste0(path_out, "/ASV_id_map.tsv"))
## Workspace
paste0(path_out, "/saved.Rdata") %>% save.image()
## Session info
sessionInfo() %>% 
  capture.output() %>% 
  writeLines(paste0(path_out, "/session.info"))

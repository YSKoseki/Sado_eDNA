# 08_Phyloseq.R
# R script to integrate ASV table, sequence data, taxonomy, phylogenetic tree,  
# sample metadata into a phyloseq object. Also computes descriptive statistics 
# of the data
# R 4.6.0

# Packages required
library(readxl); packageVersion("readxl") # 1.4.5
library(phyloseq); packageVersion("phyloseq") # 1.56.0
library(speedyseq); packageVersion("speedyseq") # 0.5.3.9021
library(tidyverse); packageVersion("tidyverse") # 2.0.0

# Paths
# Input data
path_in <- list(
  ps = "./05_Mergeruns/glob_ps.rds",
  tax = "./07_Tree/tax_tbl.rds",
  tree = "./07_Tree/treeGTRbs_tre.rds",
  smpl = "08_sampledata.xlsx",
  exclude = "08_exclude_taxa.tsv"
)
## Output directory
path_out <- "08_Phyloseq"

# Load data
ps <- path_in[["ps"]] %>% readRDS()
tax <- path_in[["tax"]] %>% readRDS()
tree <- path_in[["tree"]] %>% readRDS() %>% tidytree::as.phylo()
smpl <- path_in[["smpl"]] %>% readxl::read_excel()
exclude <- path_in[["exclude"]] %>% read_tsv()

# Convert variable types in the `smpl` dataset
smpl <- smpl %>% 
  mutate(
    type = as.factor(type), 
    year = as.factor(year), 
    siteid = as.factor(siteid), 
    riverid = as.factor(riverid), 
    coast = as.factor(coast), 
    area1 = as.factor(area1), 
    area2 = as.factor(area2),
    bottom = as.factor(bottom)
  )

# Integrate all datasets
glb_ps <- phyloseq::merge_phyloseq(
  ps,
  speedyseq::tax_table(tax),
  phyloseq::phy_tree(tree),
  speedyseq::sample_data(smpl)
)

# Human ASVs, identified by the genera listed as "human" in `exclude`
human_asv <- glb_ps %>%
  phyloseq::subset_taxa(
    genus %in% exclude$genus[exclude$reason == "human"]
  ) %>%
  taxa_names()

# Summary of read counts by sample type
reads_by_smpltype <- glb_ps %>% 
  speedyseq::psmelt(as = "tbl_df") %>% 
  select(type, Sample, Abundance) %>% 
  rename(Type = type) %>%
  summarize(
    Type = unique(Type),
    Abundance = sum(Abundance),
    .by = Sample
  ) %>% 
  summarize(
    N_samples = n(),
    Mean_reads = mean(Abundance),
    Min_reads = min(Abundance),
    Max_reads = max(Abundance),
    Total_reads = sum(Abundance),
    .by = Type
  ) %>% 
  arrange(desc(Total_reads))
dir.create(path_out, recursive = TRUE)
reads_by_smpltype %>% 
  mutate(Mean_reads = format(Mean_reads, scientific = 999, digits = 2)) %>% 
  write_tsv(paste0(path_out, "/reads_by_smpltype.tsv"))

# Examine taxa with high read counts in the negative controls of each year
# (FNC: field negative controls, LNC: laboratory negative controls)
summarize_reads_by_asv_in_ctrl <- function(ps, ctrl_type) {
  ps %>%
    speedyseq::psmelt(as = "tbl_df") %>%
    filter(type == ctrl_type) %>%
    select(type, year, OTU, Sample, nam1, Abundance) %>%
    rename(Type = type, Year = year, ASV = OTU, Taxonomy = nam1) %>%
    summarize(
      Type = unique(Type),
      Year = unique(Year),
      Taxonomy = unique(Taxonomy),
      Abundance = unique(Abundance),
      .by = c(Sample, ASV)
    ) %>%
    summarize(
      Type = unique(Type),
      Taxonomy = unique(Taxonomy),
      N_samples = n(),
      N_occur = sum(Abundance > 0),
      Mean_reads = mean(Abundance),
      min_reads = min(Abundance),
      Min_detect_reads = ifelse(
        any(Abundance > 0), min(Abundance[Abundance > 0]), NA
      ),
      Max_reads = max(Abundance),
      Total_reads = sum(Abundance),
      .by = c(Year, ASV)
    ) %>%
    mutate(
      Percent_occur = 100 * N_occur / N_samples,
      Rel_abund = 100 * Total_reads / sum(Total_reads),
      .by = Year
    ) %>%
    filter(Total_reads > 0) %>%
    select(
      Type, Year, ASV, Taxonomy, N_samples, N_occur, Percent_occur,
      Mean_reads, min_reads, Min_detect_reads, Max_reads, Total_reads, Rel_abund
    ) %>%
    arrange(Year, desc(Max_reads), ASV)
}
reads_by_asv_in_fnc <- glb_ps %>% summarize_reads_by_asv_in_ctrl("FNC")
reads_by_asv_in_lnc <- glb_ps %>% summarize_reads_by_asv_in_ctrl("LNC")
list(FNC = reads_by_asv_in_fnc, LNC = reads_by_asv_in_lnc) %>%
  iwalk(\(x, ctrl_type) {
    x %>%
      mutate(
        Percent_occur = format(Percent_occur, scientific = 999, digits = 2),
        Mean_reads = format(Mean_reads, scientific = 999, digits = 2),
        Rel_abund = format(Rel_abund, scientific = 999, digits = 2)
      ) %>%
      write_tsv(paste0(path_out, "/reads_by_ASV_in_", ctrl_type, ".tsv"))
  })

# Reads to subtract from each sample for each year and ASV: the maximum reads
# per sample in FNC plus that in LNC of the same year
reads_to_subtract <- full_join(
  reads_by_asv_in_fnc %>%
    select(Year, ASV, Taxonomy, Max_reads_FNC = Max_reads),
  reads_by_asv_in_lnc %>%
    select(Year, ASV, Taxonomy, Max_reads_LNC = Max_reads),
  by = c("Year", "ASV", "Taxonomy")
) %>%
  replace_na(list(Max_reads_FNC = 0, Max_reads_LNC = 0)) %>%
  mutate(Subtract = Max_reads_FNC + Max_reads_LNC) %>%
  arrange(Year, desc(Subtract), ASV)
reads_to_subtract %>% write_tsv(paste0(path_out, "/reads_to_subtract.tsv"))

# Subtract the reads from all samples of the same year; negative values become 0
otu <- glb_ps %>% speedyseq::otu_table() %>% as("matrix")
if (!taxa_are_rows(glb_ps)) otu <- t(otu)
smpl_year <- speedyseq::sample_data(glb_ps)$year %>% 
  as.character() %>% 
  setNames(sample_names(glb_ps))
smpl_year <- smpl_year[colnames(otu)]
for (i in seq_len(nrow(reads_to_subtract))) {
  r <- reads_to_subtract[i, ]
  in_year <- smpl_year == as.character(r$Year)
  otu[r$ASV, in_year] <- otu[r$ASV, in_year] - r$Subtract
}
otu[otu < 0] <- 0
glb_ps_sub <- glb_ps
otu_table(glb_ps_sub) <- otu_table(otu, taxa_are_rows = TRUE)

# Update the summary of read counts by sample type, excluding the human ASV
reads_by_smpltype_nh <- glb_ps %>% 
  speedyseq::psmelt(as = "tbl_df") %>% 
  filter(!OTU %in% human_asv) %>%
  select(type, Sample, Abundance) %>% 
  rename(Type = type) %>%
  summarize(
    Type = unique(Type),
    Abundance = sum(Abundance),
    .by = Sample
  ) %>% 
  summarize(
    N_samples = n(),
    Mean_reads = mean(Abundance),
    Min_reads = min(Abundance),
    Min_detect_reads = ifelse(
      any(Abundance > 0), min(Abundance[Abundance > 0]), NA
    ),
    Max_reads = max(Abundance),
    Total_reads = sum(Abundance),
    .by = Type
  ) %>% 
  arrange(desc(Total_reads))
reads_by_smpltype_nh %>% 
  mutate(Mean_reads = format(Mean_reads, scientific = 999, digits = 2)) %>% 
  write_tsv(paste0(path_out, "/reads_by_smpltype_nh.tsv"))

# Summary of read counts for each ASV in field samples
reads_by_asv_in_smpl <- glb_ps %>% 
  phyloseq::subset_samples(type == "SMPL") %>% 
  speedyseq::psmelt(as = "tbl_df") %>% 
  select(OTU, Abundance, nam1) %>% 
  rename(ASV = OTU, Taxonomy = nam1) %>% 
  #mutate(Occurrence = if_else(Abundance > 0, 1L, 0L)) %>% 
  summarize(
    Taxonomy = unique(Taxonomy),
    N_samples = n(),
    N_occur = sum(Abundance > 0),
    Mean_reads = mean(Abundance),
    Min_reads = min(Abundance),
    Min_detect_reads = ifelse(
      any(Abundance > 0), min(Abundance[Abundance > 0]), NA
    ),
    Max_reads = max(Abundance),
    Total_reads = sum(Abundance),
    .by = ASV
  ) %>% 
  mutate(
    Percent_occur = 100 * N_occur / N_samples,
    Rel_abund = 100 * Total_reads / sum(Total_reads)
  ) %>% 
  select(
    ASV, Taxonomy, N_samples, N_occur, Percent_occur, 
    Mean_reads, Min_reads, Min_detect_reads, Max_reads, Total_reads, Rel_abund
  ) %>% 
  arrange(ASV)
reads_by_asv_in_smpl %>% 
  mutate(
    Percent_occur = format(Percent_occur, scientific = 999, digits = 2),
    Mean_reads = format(Mean_reads, scientific = 999, digits = 2),
    Rel_abund = format(Rel_abund, scientific = 999, digits = 2)
  ) %>% 
  write_tsv(paste0(path_out, "/reads_by_ASV_in_SMPL.tsv"))

# Create a phyloseq object for freshwater fish
## Use the reads after subtracting the negative-control reads, and exclude
## human and marine taxa by genus, as listed in `exclude`
fw_ps <- glb_ps_sub %>%
  phyloseq::subset_taxa(!genus %in% exclude$genus)
## ASVs excluded
excluded_asv <- glb_ps %>%
  speedyseq::tax_table() %>%
  as.data.frame() %>%
  rownames_to_column("ASV") %>%
  filter(!ASV %in% taxa_names(fw_ps)) %>%
  select(ASV, genus, nam1) %>%
  left_join(exclude, by = "genus")
excluded_asv %>% write_tsv(paste0(path_out, "/excluded_ASV.tsv"))

# Subset the phyloseq object to field samples only
fw_smpl_ps <- fw_ps %>% 
  phyloseq::subset_samples(type == "SMPL")
　
# Outputs
## R objects
saveRDS(fw_ps, paste0(path_out, "/fw_ps.rds"))
saveRDS(fw_smpl_ps, paste0(path_out, "/fw_smpl_ps.rds"))
## Workspace
paste0(path_out, "/saved.Rdata") %>% save.image()
## Session info
sessionInfo() %>% 
  capture.output() %>% 
  writeLines(paste0(path_out, "/session.info"))

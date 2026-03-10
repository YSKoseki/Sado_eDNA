# 08_Phyloseq.R
# R script to integrate ASV table, sequence data, taxonomy, phylogenetic tree,  
# sample metadata into a phyloseq object. Also computes descriptive statistics 
# of the data
# R 4.5.2

# Packages required
library(readxl); packageVersion("readxl") # 1.4.5
library(phyloseq); packageVersion("phyloseq") # 1.54.2
library(speedyseq); packageVersion("speedyseq") # 0.5.3.9021
library(tidyverse); packageVersion("tidyverse") # 2.0.0

# Paths
# Input data
path_in <- list(
  ps = "./05_Mergeruns/glob_ps.rds",
  tax = "./07_Tree/tax_tbl.rds",
  tree = "./07_Tree/treeGTRbs_tre.rds",
  smpl = "08_sampledata.xlsx"
)
## Output directory
path_out <- "08_Phyloseq"

# Load data
ps <- path_in[["ps"]] %>% readRDS()
tax <- path_in[["tax"]] %>% readRDS()
tree <- path_in[["tree"]] %>% readRDS() %>% tidytree::as.phylo()
smpl <- path_in[["smpl"]] %>% readxl::read_excel()

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

# Examine taxa with high read counts in the field negative controls
reads_by_asv_in_fnc <- glb_ps %>% 
  speedyseq::psmelt(as = "tbl_df") %>% 
  filter(type == "FNC") %>% 
  select(type, OTU, Sample, nam1, Abundance) %>% 
  rename(Type = type, ASV = OTU, Taxonomy = nam1) %>%
  summarize(
    Type = unique(Type),
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
    .by = ASV
  ) %>%
  mutate(
    Percent_occur = 100 * N_occur / N_samples,
    Rel_abund = 100 * Total_reads / sum(Total_reads)
  ) %>% 
  filter(Total_reads > 0) %>% 
  select(
    Type, ASV, Taxonomy, N_samples, N_occur, Percent_occur, 
    Mean_reads, min_reads, Min_detect_reads, Max_reads, Total_reads, Rel_abund
  ) %>% 
  arrange(desc(Max_reads), ASV)
reads_by_asv_in_fnc %>% 
  mutate(
    Percent_occur = format(Percent_occur, scientific = 999, digits = 2),
    Mean_reads = format(Mean_reads, scientific = 999, digits = 2),
    Rel_abund = format(Rel_abund, scientific = 999, digits = 2)
  ) %>% 
  write_tsv(paste0(path_out, "/reads_by_ASV_in_FNC.tsv"))

# Use the maximum non-human read count as the ASV filtering cutoff
max_reads_in_fnc <- reads_by_asv_in_fnc %>% 
  filter(Taxonomy != "Homo sapiens") %>%
  select(Max_reads) %>% 
  pull() %>% 
  max()

# Update the summary of read counts by sample type, excluding the human ASV
reads_by_smpltype_nh <- glb_ps %>% 
  speedyseq::psmelt(as = "tbl_df") %>% 
  filter(nam1 != "Homo sapiens") %>% 
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
human <- c("ASV069")
names(human) <- c("Homo sapiens")
marine <- c("ASV032", "ASV045", "ASV057", "ASV068", "ASV074", "ASV096")
names(marine) <- c(
  "Seriola dumerili", "Seriola quinqueradiata", "Epinephelus akaara", 
  "Seriola quinqueradiata", "Epinephelus akaara", "Sebastes sp."
)
fw_ps <- glb_ps %>% 
  taxa_names() %>% 
  `[`(!. %in% c(human, marine)) %>% 
  phyloseq::prune_taxa(glb_ps)

# Subset the phyloseq object to field samples only
fw_smpl_ps <- fw_ps %>% 
  phyloseq::subset_samples(type == "SMPL")
　
# Outputs
## R objects
saveRDS(fw_ps, paste0(path_out, "/fw_ps.rds"))
saveRDS(fw_smpl_ps, paste0(path_out, "/fw_smpl_ps.rds"))
saveRDS(max_reads_in_fnc, paste0(path_out, "/max_reads_in_fnc.rds"))
## Workspace
paste0(path_out, "/saved.Rdata") %>% save.image()
## Session info
sessionInfo() %>% 
  capture.output() %>% 
  writeLines(paste0(path_out, "/session.info"))

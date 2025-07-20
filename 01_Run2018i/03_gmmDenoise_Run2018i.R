# 03_gmmDenoise_Run2018i.R
# R script to filter false-positive ASVs using gmmDenoise (Koseki et al. 2025)
# Reference https://github.com/YSKoseki/gmmDenoise
# R 4.4.2

# Packages required
library(gmmDenoise); packageVersion("gmmDenoise") # 0.3.1
library(cowplot); packageVersion("cowplot") # 1.1.3
library(scales); packageVersion("scales") # 1.4.0
library(tidyverse); packageVersion("tidyverse") # 2.0.0

# Paths for:
## ASV table
path_asvtab <- "02_DADA2_Run2018i/05_saved_rds/ASV_tab.rds"
## Summary table of the preceding Cutadapt–DADA2 pipeline
path_sumtab <- "02_DADA2_Run2018i/05_saved_rds/summary_tab.rds"
## Output directory
path_out <- "03_gmmDenoise_Run2018i"

# Load data
asv_tab <- readRDS(path_asvtab)
summary_tab <- readRDS(path_sumtab)

# Create ASV read count vector
read_count <- asv_tab %>% colSums()

# gmmDenoise filtering
## Step1: Check ASV read count distribution (log-scale histogram)
(fig_hist <- read_count %>% asvhist())
## Step 2: Cross-validation for selecting the number of mixture components
set.seed(101)
crossval <- read_count %>% log10() %>% gmmcv(maxit = 5000)
maxlik_id <- crossval$log.lik %>% which.max()
k <- crossval$k[maxlik_id]
(fig_cv <- crossval %>% autoplot())
# Step 3: Model fitting with the selected number of mixture components
set.seed(101)
gmm_fit <- read_count %>% log10() %>% gmmem(k, maxit = 5000)
# Step 4: Determine cut-off threshold (upper one-sided 95% CL of the second
#  uppermost component)
ln_thresh <- gmm_fit %>% quantile()
thresh <- ceiling(10^ln_thresh)
fig_gmm <- gmm_fit %>% 
  autoplot(vline = ln_thresh) + 
  annotate(
    "text",
    label = paste0("Threshold = ", round(ln_thresh, 2)),
    size = 5,
    x = Inf, y = Inf,
    hjust = 1.3, vjust = 2
  )
# Step 5: Plot the gmmDenoise results
dir.create(path_out, recursive=TRUE)
theme_set(cowplot::theme_cowplot())
figs <- plot_grid(
  fig_hist, fig_cv, fig_gmm,
  labels = c("(a)", "(b)", "(c)"), 
  nrow = 3,
  align = "hv",
  label_x = .13, label_y = .97
)
save_plot(
  paste0(path_out, "/Fig_gmm.svg"),
  figs,
  ncol = 1,
  nrow = 3
)
# Step 6: Filter ASVs with the cut-off threshold
retained <- which(read_count > thresh)
read_count_flt <- read_count[retained]
asv_tab_flt <- asv_tab[, retained]

# Create a summary table for the Cutadapt–DADA2–gmmDenoise pipeline
read_count_tab <- tibble(
  sample = asv_tab_flt %>% rowSums() %>% names(),
  gmmdenoise = asv_tab_flt %>% rowSums()
)
summary_tab_flt <- summary_tab %>% 
  select(-pct_remain) %>% 
  left_join(read_count_tab) %>% 
  mutate(gmmdenoise = replace_na(gmmdenoise, 0)) %>% 
  mutate(pct_remain = percent(gmmdenoise / initial, accuracy = 0.1))
summary_tab_flt %>% 
  print(n = summary_tab_flt %>% nrow())
write_tsv(summary_tab_flt, file = paste0(path_out, "/gmmDenoise_sum.tsv"))

# Outputs
## R objects
saveRDS(asv_tab_flt, paste0(path_out, "/ASV_tab_gmm.rds"))
## Workspace
paste0(path_out, "/saved.Rdata") %>% save.image()
## Session info
sessionInfo() %>% 
  capture.output() %>% 
  writeLines(paste0(path_out, "/session.info"))


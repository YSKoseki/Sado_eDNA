# 07_Tree.R
# R script to construct and plot phylogenetic trees
# R 4.4.2

# Packages required
library(data.table); packageVersion("data.table") # 1.17.0
library(Biostrings); packageVersion("Biostrings") # 2.74.1
library(DECIPHER); packageVersion("DECIPHER") # 3.2.0
library(phangorn); packageVersion("phangorn") # 2.12.1
library(parallel); packageVersion("parallel") # 4.4.2
library(tidytree); packageVersion("tidytree") # 0.4.6
library(ggtree); packageVersion("ggtree") # 3.14.0
library(cowplot); packageVersion("cowplot") # 1.1.3
library(tidyverse); packageVersion("tidyverse") # 2.0.0

# Paths
# Input data
path_in <- list(
  ps = "./05_Mergelib/glob_ps.rds",
  cl = "./06_Claident/taxonomy_merger.tsv"
)
## Output directory
path_out <- "07_Tree"

# Load data
ps <- path_in$ps %>% readRDS()
cl <- path_in$cl %>% fread() %>% as_tibble()

# Create ASV names from Claident-assigned taxonomy
tax_tbl <- cl %>% 
  select(query, class, order, family, genus, species) %>% 
  rename(label = query) %>% 
  mutate(
    nam1 = if_else(
      species != "", species,
      if_else(
        genus != "", paste0(genus, " sp."), 
        if_else(
          family != "", paste0(family, " sp."), 
          if_else(
            order != "", paste0(order, " sp."), 
            if_else(
              class != "", paste0(class, " sp."),
              "Unassigned"
            )
          )
        )
      )
    )
  ) %>% 
  mutate(nam2 = paste0(nam1, " [", label, "]")) %>% 
  mutate(sp_assign = if_else(species == "", "no", "yes")) %>% 
  mutate(seq = refseq(ps) %>% as.character())

# Construct phylogenetic trees
#  See `vignette("Trees", package = "phangorn")`
## Align sequences
alignment <- refseq(ps) %>% 
  DECIPHER::AlignSeqs(anchor = NA, processors = NULL, verbose = FALSE) 
## Convert the alignment object into a phyDat object
align_pha <- alignment %>% as.matrix() %>% phangorn::phyDat(type = "DNA")
## Compute pairwise distances of aligned sequences
align_dis <- align_pha %>% phangorn::dist.ml()
## Build a tree using the neighbor-joining method
tree_NJ <- align_dis %>% phangorn::NJ()
## Build a tree using the maximum likelihood method with the GTR+I+G model (the 
##  General Time Reversible model with corrections for invariant sites and 
##  gamma-distributed rate variation)
tree_GTR <- tree_NJ %>% 
  phangorn::pml(data = align_pha) %>% 
  update(k = 4, inv = 0.2) %>%
  phangorn::optim.pml(
    model = "GTR",
    optInv = TRUE,
    optGamma = TRUE,
    rearrangement = "stochastic",
    control = phangorn::pml.control(trace = 0)
  )
## Compute bootstrap support values for the ML tree
set.seed(123)
bs <- tree_GTR %>% 
  phangorn::bootstrap.pml(
    bs = 1000, 
    model = "GTR",
    optInv = TRUE,
    optGamma = TRUE,
    rearrangement = "stochastic",
    multicore = TRUE,
    mc.cores = parallel::detectCores()
  )
## The ML tree with bootstrap support values
tree_GTRbs <- tree_GTR
tree_GTRbs$tree <- tree_GTR$tree %>% phangorn::plotBS(bs)

# Add detailed annotation data to the ML tree
#  Reference: https://yulab-smu.top/treedata-book/index.html
## Merge taxonomic information into the above tree data
treeGTRbs_tbl <- tree_GTRbs$tree %>% 
  as_tibble() %>% 
  mutate(
    bsval = case_when(
      !str_detect(label, "^ASV") ~ label, 
      TRUE ~ NA_character_
    ),
    bsval = as.numeric(bsval),
  ) %>% 
  full_join(tax_tbl, by = "label")
## Convert the tibble object into treedata object
treeGTRbs_tre <- treeGTRbs_tbl %>% tidytree::as.treedata()
## Remove the root node 69 (Homo sapiens [ASV069])
treeGTRbs_tre2 <- treeGTRbs_tre %>% tidytree::drop.tip(69)

# Plot the global tree
dir.create(path_out, recursive = TRUE)
param <- list(
  linesize = .2, textsize = 1.8, nodesize = 1.2, bsvsize = 1.2, anntsize = 6
)
ggt_glb <- treeGTRbs_tre2 %>% 
  ggtree::ggtree(size = param$linesize) +
  # Annotate tree tips with ASV names
  ggtree::geom_text2(
    aes(subset = isTip, label = nam2), size = param$textsize, hjust = -.08
  ) +
  # Annotate internal nodes with their ID numbers
  ggtree::geom_text2(
    aes(subset = !isTip, label = node), 
    size = param$nodesize, 
    colour = "blue",
    hjust = -.1
  ) +
  # Annotate branches with bootstrap support values if they are >0.7
  geom_text2(
    aes(
      subset = !isTip,
      x = branch,
      label = if_else(bsval > 0.7, sprintf("%.2f", bsval), NA_character_),
      fontface = 3
    ),
    size = param$bsvsize,
    vjust = -.5
  ) +
  scale_x_continuous(limits = c(0, 1.5)) +
  ggtree::geom_treescale(
    x = 0, y = -1.5, 
    width = .1,
    linesize = param$linesize, 
    fontsize = param$textsize
  ) +
  annotate(
    "text", 
    x = 0.3, y = 102, 
    label = "Global tree with all fish ASVs", 
    size = param$anntsize
  )
ggt_glb
## Save the plot
ggsave(
  paste0(path_out, "/Globaltree.svg"), 
  device = "svg", bg = "white",
  width = 21, height = 25, units = "cm"
)

# Plot subtrees for inspection
## Custom function to plot a subtree
subtree <- function(subsetnode, collapsenode = NULL, levels_back = 0, xlim = c(0, 1), linesize = .1, textsize = 1, scalewidth = .01, title = "", titlesuffix = "clade", titlesize = 3.5, width = 4, height = 3) {
  p <- treeGTRbs_tre2 %>% 
    tidytree::tree_subset(node = subsetnode, levels_back = levels_back) %>% 
    ggtree::ggtree(size = linesize) +
    # Annotate tree tips with ASV names
    ggtree::geom_text2(
      aes(subset = isTip, label = nam2), 
      size = textsize, 
      hjust = -.08
    ) +
    scale_x_continuous(limits = xlim) +
    ggtree::geom_treescale(
      x = 0, y = -.5, 
      width = scalewidth, 
      linesize = linesize, 
      fontsize = textsize
    ) +
    labs(title = paste0("The ", title, titlesuffix)) +
    theme(text = element_text(size = titlesize))
  if (is.null(collapsenode)) {
    p
  } else {
    p2 <- p %>% 
      collapse(node = collapsenode) +
      geom_point2(
        aes(subset = (node == collapsenode)), 
        shape = 18, size = .3, fill = "black"
      )
    p2
  }
  ggsave(
    paste0(path_out, paste0("/", title, ".svg")), 
    device = "svg", bg = "white",
    width = width, height = height, units = "cm"
  )
}
subtree(subsetnode = 151, xlim = c(0, .25), title = "Misgurnus")
subtree(subsetnode = 162, xlim = c(0, .03), title = "Carassius-Cyprinus")
subtree(subsetnode = 189, xlim = c(0, .15), title = "Cottus")
subtree(subsetnode = 195, xlim = c(0, .22), title = "Gymnogobius-Luciogobius")
subtree(
  subsetnode = 106, collapsenode = 129,
  xlim = c(0, .077), textsize = 0.7,
  title = "Rhinogobius", titlesuffix = " branches"
)

# Outputs
## R objects
saveRDS(tax_tbl, paste0(path_out, "/tax_tbl.rds"))
saveRDS(treeGTRbs_tre, paste0(path_out, "/treeGTRbs_tre.rds"))
## Workspace
paste0(path_out, "/saved.Rdata") %>% save.image()
## Session info
sessionInfo() %>% 
  capture.output() %>% 
  writeLines(paste0(path_out, "/session.info"))

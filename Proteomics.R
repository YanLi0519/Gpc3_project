# ==============================================================================
# Proteomics analysis pipeline (SAT vs VAT; Lean vs Obese)
# Author: Yan Li
# Contact: liyan@stu.pku.edu.cn
# Lab: CarlosLab@PKU
# License: MIT (recommended for GitHub) or CC BY 4.0 for code+docs
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggVennDiagram)
  library(VennDiagram)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ReactomePA)
  library(pheatmap)
  library(here)
})

source("R/00_utils.R")
source("R/01_load_clean.R")
source("R/02_density.R")
source("R/03_enriched_sets.R")
source("R/04_venn.R")
source("R/05_lollipop.R")
source("R/06_enrichment_GO.R")
source("R/07_obese_specific.R")

# --------------------------
# Parameters
# --------------------------
params <- list(
  input_csv = here("data", "all_sample.csv"),
  out_fig   = here("results", "figures"),
  out_tbl   = here("results", "tables"),
  log2_pseudo = 1,          # log2(x + 1)
  fc_threshold = 0.3,
  density_xlim = c(-1.5, 1.5),
  density_ylim = c(0, 3),
  top_n_lollipop = 20,
  go_ont = "BP",
  go_show = 20
)

dir_create(params$out_fig)
dir_create(params$out_tbl)

# --------------------------
# 1) Load & clean
# --------------------------
dat <- load_and_clean_proteomics(
  file = params$input_csv,
  log2_pseudo = params$log2_pseudo,
  dedup_method = "least_na"  # "first" / "mean" / "max" / "least_na"
)

all_sample <- dat$df
meta <- dat$meta

# --------------------------
# 2) Density plots per individual
# --------------------------
p_density <- plot_density_all_individuals(
  df = all_sample,
  sample_ids = meta$sample_ids,
  xlim = params$density_xlim,
  ylim = params$density_ylim
)

save_ggplot_pdf(
  plot = p_density,
  filename = file.path(params$out_fig, "SAT_VAT_log2FC_density_per_individual.pdf"),
  width = 14, height = 7
)

# --------------------------
# 3) Enriched sets (per person) + intersection
# --------------------------
sets <- get_enriched_sets_grouped(
  df = all_sample,
  lean_ids = meta$lean_ids,
  obese_ids = meta$obese_ids,
  threshold = params$fc_threshold
)

# save intersections
write_csv(tibble(Gene = sets$lean$SAT_intersect), file.path(params$out_tbl, "lean_SAT_enriched.csv"))
write_csv(tibble(Gene = sets$lean$VAT_intersect), file.path(params$out_tbl, "lean_VAT_enriched.csv"))
write_csv(tibble(Gene = sets$obese$SAT_intersect), file.path(params$out_tbl, "obese_SAT_enriched.csv"))
write_csv(tibble(Gene = sets$obese$VAT_intersect), file.path(params$out_tbl, "obese_VAT_enriched.csv"))

# --------------------------
# 4) Venn plots (center-only or full)
# --------------------------
p_venn_center <- plot_venn_4panels_center_only(sets)
save_ggplot_pdf(p_venn_center, file.path(params$out_fig, "Venn_center_only_4panels.pdf"), 10, 8)

# --------------------------
# 5) Lollipop (Top N by p-value)
# --------------------------
p_lollipop_lean  <- plot_lollipop_topN(all_sample, group_prefix = "L",
                                       gene_lists = list(SAT=sets$lean$SAT_intersect, VAT=sets$lean$VAT_intersect),
                                       top_n = params$top_n_lollipop)

p_lollipop_obese <- plot_lollipop_topN(all_sample, group_prefix = "O",
                                       gene_lists = list(SAT=sets$obese$SAT_intersect, VAT=sets$obese$VAT_intersect),
                                       top_n = params$top_n_lollipop)

save_ggplot_pdf(p_lollipop_lean,  file.path(params$out_fig, "Lean_Top20_Lollipop.pdf"), 14, 3)
save_ggplot_pdf(p_lollipop_obese, file.path(params$out_fig, "Obese_Top20_Lollipop.pdf"), 14, 3)

# --------------------------
# 6) GO enrichment (BP)
# --------------------------
run_go_enrichment_all(
  sets = sets,
  out_dir = params$out_tbl,
  fig_dir = params$out_fig,
  ont = params$go_ont,
  showCategory = params$go_show
)

# --------------------------
# 7) Obese-specific + heatmaps (save properly)
# --------------------------
run_obese_specific_heatmaps(
  df = all_sample,
  meta = meta,
  sets = sets,
  out_dir = params$out_fig
)

# --------------------------
# Record session info for reproducibility
# --------------------------
writeLines(capture.output(sessionInfo()), file.path(params$out_tbl, "sessionInfo.txt"))
message("Done. Outputs in: ", here("results"))
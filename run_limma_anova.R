# R/run_limma_anova.R
# Moderated ANOVA-like testing across groups using limma (one model for all metabolites).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(limma)
  library(ggplot2)
})

dat <- read_csv("data/demo_data.csv", show_col_types = FALSE)

dat <- dat %>% mutate(Group = factor(Group, levels = c("Group1","Group2","Group3","Group4")))

meta_cols <- setdiff(names(dat), c("sample_ID","Group"))
M <- dat %>% select(all_of(meta_cols)) %>% as.matrix() %>% t()

for (i in seq_len(nrow(M))) {
  if (anyNA(M[i, ])) {
    mu <- mean(M[i, ], na.rm = TRUE)
    M[i, is.na(M[i, ])] <- mu
  }
}

keep <- apply(M, 1, function(x) sd(x) > 0)
M <- M[keep, , drop = FALSE]
features <- rownames(M)

grp <- dat$Group
design <- model.matrix(~ 0 + grp)
colnames(design) <- levels(grp)

fit <- lmFit(M, design)
contrast_mat <- makeContrasts(
  G2_vs_G1 = Group2 - Group1,
  G3_vs_G1 = Group3 - Group1,
  G4_vs_G1 = Group4 - Group1,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

extract_contrast <- function(fit2, coef_name, n = Inf) {
  tt <- topTable(fit2, coef = coef_name, number = n, sort.by = "P")
  tt %>% rownames_to_column("Metabolite") %>%
    mutate(Contrast = coef_name) %>%
    select(Contrast, Metabolite, logFC, AveExpr, t, P.Value, adj.P.Val, B)
}

res <- bind_rows(
  extract_contrast(fit2, "G2_vs_G1"),
  extract_contrast(fit2, "G3_vs_G1"),
  extract_contrast(fit2, "G4_vs_G1")
)

dir.create("results", showWarnings = FALSE)
readr::write_csv(res, "results/limma_anova_results.csv")
message("Wrote results/limma_anova_results.csv")

dir.create("figs", showWarnings = FALSE)
top_hit <- res %>% filter(Contrast == "G2_vs_G1") %>% arrange(adj.P.Val) %>% slice_head(n = 1) %>% pull(Metabolite)

if (length(top_hit) == 1) {
  p <- ggplot(dat, aes(x = Group, y = .data[[top_hit]])) +
    geom_boxplot() +
    labs(title = paste("Top metabolite (G2 vs G1):", top_hit), y = top_hit) +
    theme_minimal()
  ggsave(filename = file.path("figs", paste0("boxplot_", top_hit, ".png")),
         plot = p, width = 6, height = 4, dpi = 150)
  message("Saved figs/boxplot_", top_hit, ".png")
}

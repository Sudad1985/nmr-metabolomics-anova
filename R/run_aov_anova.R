# R/run_aov_anova.R
# Base R ANOVA per metabolite with BH correction across metabolites.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(broom)
})

dat <- read_csv("data/demo_data.csv", show_col_types = FALSE) %>%
  mutate(Group = factor(Group, levels = c("Group1","Group2","Group3","Group4")))

meta_cols <- setdiff(names(dat), c("sample_ID","Group"))

do_aov <- function(y, g) {
  df <- data.frame(y = y, g = g)
  df <- df[is.finite(df$y), ]
  if (length(unique(df$g)) < 2) return(NULL)
  m <- aov(y ~ g, data = df)
  broom::tidy(m) %>% filter(term == "g")
}

aov_res <- lapply(meta_cols, function(met) {
  out <- do_aov(dat[[met]], dat$Group)
  if (is.null(out)) return(NULL)
  out$Metabolite <- met
  out
}) %>%
  bind_rows() %>%
  relocate(Metabolite)

aov_res <- aov_res %>%
  mutate(adj.p.value = p.adjust(p.value, method = "BH"))

dir.create("results", showWarnings = FALSE)
readr::write_csv(aov_res, "results/aov_anova_results.csv")
message("Wrote results/aov_anova_results.csv")

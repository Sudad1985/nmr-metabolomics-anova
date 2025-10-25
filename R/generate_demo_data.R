# R/generate_demo_data.R
# Creates a synthetic metabolomics table: samples in rows, metabolites in columns, plus a Group factor.

set.seed(42)

n_per_group <- 15
groups <- factor(rep(paste0("Group", 1:4), each = n_per_group),
                 levels = paste0("Group", 1:4))

# A small panel of metabolites
metabolites <- c(
  "X2_Aminobutyric_acid","Alanine","Asparagine","Creatinine","Glutamic_acid",
  "Glutamine","Glycine","Histidine","Isoleucine","Leucine","Lysine",
  "Phenylalanine","Tyrosine","Valine","X2_Hydroxybutyric_acid","Acetic_acid",
  "Citric_acid","Lactic_acid","X3_Hydroxybutyric_acid","Acetoacetic_acid",
  "Acetone","Pyruvic_acid","Glucose"
)

n <- length(groups)
dat <- data.frame(
  sample_ID = paste0("S", seq_len(n)),
  Group = groups,
  stringsAsFactors = FALSE
)

for (met in metabolites) {
  base_mean <- runif(1, 5, 12)
  group_shift <- c(0, rnorm(1, 0.8, 0.4), rnorm(1, -0.6, 0.4), rnorm(1, 0.4, 0.4))
  vals <- numeric(n)
  for (g in 1:4) {
    idx <- which(groups == paste0("Group", g))
    vals[idx] <- rnorm(length(idx), mean = base_mean + group_shift[g], sd = runif(1, 0.5, 1.2))
  }
  dat[[met]] <- round(vals, 3)
}

dir.create("data", showWarnings = FALSE)
write.csv(dat, "data/demo_data.csv", row.names = FALSE)
message("Wrote data/demo_data.csv")

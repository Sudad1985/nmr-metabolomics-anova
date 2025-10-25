# NMR Metabolomics — Group ANOVA in R (limma & base R)

Reproducible R scripts for comparing metabolite levels across groups using:
- **limma** (moderated ANOVA style with contrasts | recommended)
- **base R `aov()`** with BH multiple-testing correction

> This code generalizes analysis I performed during a 2024 research internship
> at the Institute of Functional Genomics, University of Regensburg.  
> **Note:** The data here are synthetic demo data; original research data are not shared.

## How to run
```bash
Rscript R/generate_demo_data.R       # writes data/demo_data.csv
Rscript R/run_limma_anova.R          # writes results/limma_anova_results.csv
Rscript R/run_aov_anova.R            # writes results/aov_anova_results.csv
```

## Outputs
- `results/limma_anova_results.csv` — limma results for G2/G3/G4 vs G1
- `results/aov_anova_results.csv` — base-R ANOVA results (one row per metabolite)
- `figs/boxplot_<met>.png` — example plot for the top significant metabolite (optional)

## Stack
R, limma, dplyr, broom, ggplot2

## Related publication (preprint)
If your internship supported a medRxiv preprint, cite it here with a DOI.  
*I contributed statistical analysis and R code but I am not listed as an author.*

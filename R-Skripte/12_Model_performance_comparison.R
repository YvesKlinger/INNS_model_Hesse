# Model_performance_comparison.R
#
# Statistical comparison of SDM test performance (AUC, TSS) between data
# sources (GBIF, HLNUG, combined) of the local models across all individual
# model runs (iterations).
#
# All models are beta regressions (glmmTMB, logit link). For each metric:
#   1. Overall effect of data source, accounting for species identity:
#        metric ~ source + species   vs.   metric ~ species
#   2. Does the effect of data source depend on the species?
#        metric ~ source * species   vs.   metric ~ source + species
#   3. Comparison of data sources within each species (10 runs per source):
#        metric ~ source   vs.   metric ~ 1
#      with pairwise differences as estimated marginal means on the response
#      scale with 95% CIs (emmeans, Tukey-adjusted).
# Effects are tested with likelihood-ratio tests. Model assumptions are
# checked with simulated quantile residuals (DHARMa) for every fitted model.
#
# Also produces the AUC/TSS boxplots of the Supplementary material (App. 4.1).
#
# Input:  per-run metric files written by 10_SDM_Code_Model_final.R
#         (Model_results/Metrics/{AUC,TSS}/<metric>__local_<origin>_<species>_<run>.csv).
#         If Model_results/Metrics/ does not exist, the archived per-run
#         metrics in "Model result metrics/Metrics.zip" are used instead.
# Output: Model_results/Model_comparison/

#Setup
##########################################################################
##########################################################################

#Packages
library(dplyr)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(flextable)
library(officer)
library(ggplot2)

# Set the working directory to the project root (the folder that contains
# R-Skripte/, Model_results/, "Model result metrics/", etc.).
# Recommended: open this project via the provided .Rproj file, which sets
# the working directory to the project root automatically. Otherwise:
# setwd("path/to/INNS_model_Hesse")

out_dir <- "Model_results/Model_comparison/"
dir.create(file.path(out_dir, "DHARMa"), recursive = TRUE, showWarnings = FALSE)

set.seed(42)


#Load per-run metrics
##########################################################################
##########################################################################

metrics_dir <- "Model_results/Metrics"
if (!dir.exists(file.path(metrics_dir, "AUC"))) {
  message("Model_results/Metrics/ not found - using archived 'Model result metrics/Metrics.zip'")
  unzip_dir <- file.path(tempdir(), "Metrics_unzipped")
  unzip("Model result metrics/Metrics.zip", exdir = unzip_dir)
  metrics_dir <- file.path(unzip_dir, "Metrics")
}

read_metric <- function(metric) {
  # per-run files of the local models only (the *_mean.csv files are averages over runs)
  pattern <- paste0("^", metric, "__local_([a-z]+)_([0-9]+)_([0-9]+)\\.csv$")
  files <- list.files(file.path(metrics_dir, toupper(metric)), pattern = pattern, full.names = TRUE)
  parts <- regmatches(basename(files), regexec(pattern, basename(files)))
  data.frame(
    metric  = toupper(metric),
    origin  = sapply(parts, `[`, 2),
    species = sapply(parts, `[`, 3),
    run     = as.integer(sapply(parts, `[`, 4)),
    # files were written with write.csv2 (single column "x"); TSS files use a
    # decimal comma, AUC files (pROC auc objects) a decimal point
    value   = sapply(files, function(f) {
      as.numeric(sub(",", ".", read.csv2(f, colClasses = "character")[[1]][1]))
    }),
    row.names = NULL
  )
}

perf <- bind_rows(read_metric("auc"), read_metric("tss")) %>%
  mutate(
    source = factor(origin,
                    levels = c("combined", "gbif", "hlnug"),
                    labels = c("COM", "GBIF", "HLNUG")),
    species = factor(species,
                     levels = c("3190653", "3034825", "2891770"),
                     labels = c("A. altissima", "H. mantegazzianum", "I. glandulifera"))
  )

# beta regression requires values strictly inside (0, 1)
if (any(perf$value <= 0 | perf$value >= 1)) {
  stop("Metric values on or outside the (0, 1) boundary - apply a transformation ",
       "(e.g. Smithson & Verkuilen 2006) before fitting beta models.")
}

write.csv(perf, paste0(out_dir, "performance_per_run.csv"), row.names = FALSE)

# descriptive statistics per source (pooled over species) and per species and source
desc_stats <- bind_rows(
  perf %>% group_by(metric, source) %>% mutate(species = "all species"),
  perf %>% group_by(metric, species, source) %>% mutate(species = as.character(species))
) %>%
  group_by(metric, species, source) %>%
  summarise(n = n(), mean = mean(value), sd = sd(value), se = sd / sqrt(n),
            min = min(value), max = max(value), .groups = "drop")
write.csv(desc_stats, paste0(out_dir, "descriptive_stats.csv"), row.names = FALSE)
print(desc_stats, n = Inf)


#Boxplots of AUC and TSS per data source and species (Supplementary material, App. 4.1)
##########################################################################
##########################################################################

plot_data <- perf %>%
  mutate(model = factor(source, levels = c("COM", "GBIF", "HLNUG"), labels = c("COM", "CS", "StAg")))

for (metric in c("AUC", "TSS")) {
  p_box <- ggplot(filter(plot_data, metric == !!metric), aes(x = model, y = value, fill = model)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ species, scales = "free_x", ncol = 3) +
    labs(title = paste(metric, "Scores by Model Type and Species"),
         x = "Model", y = paste(metric, "Score"), fill = "Model") +
    theme_minimal() +
    theme(legend.position = "top", strip.text = element_text(face = "italic"))
  ggsave(paste0(out_dir, "Supplementary_figure_", metric, "_boxplot.png"), p_box,
         width = 10, height = 6, dpi = 300, bg = "white")
}


#Helper functions
##########################################################################
##########################################################################

fit_beta <- function(formula, data) {
  glmmTMB(formula, family = beta_family(link = "logit"), data = data)
}

# likelihood-ratio test of a reduced vs. a full model
lrt_row <- function(m_reduced, m_full, metric, test, species = "all species") {
  lrt <- anova(m_reduced, m_full)
  data.frame(
    metric  = metric,
    test    = test,
    species = species,
    n_obs   = nobs(m_full),
    chisq   = lrt$Chisq[2],
    df      = lrt$`Chi Df`[2],
    p       = lrt$`Pr(>Chisq)`[2],
    AIC_reduced = lrt$AIC[1],
    AIC_full    = lrt$AIC[2]
  )
}

# DHARMa residual check: diagnostic plots + test table
dharma_check <- function(model, data, model_name) {
  sim_res <- simulateResiduals(model, n = 1000)

  png(paste0(out_dir, "DHARMa/DHARMa_", model_name, ".png"), width = 2400, height = 1200, res = 200)
  plot(sim_res)
  dev.off()

  # residuals per data source + test for homogeneity of variance between sources
  png(paste0(out_dir, "DHARMa/DHARMa_", model_name, "_by_source.png"), width = 2400, height = 1200, res = 200)
  cat_test <- testCategorical(sim_res, catPred = data$source)
  dev.off()

  uniformity <- testUniformity(sim_res, plot = FALSE)
  dispersion <- testDispersion(sim_res, plot = FALSE)
  outliers   <- testOutliers(sim_res, type = "bootstrap", plot = FALSE)

  data.frame(
    model = model_name,
    test = c("KS uniformity", "Dispersion", "Outliers", "Homogeneity of variance between sources (Levene)"),
    statistic = c(uniformity$statistic, dispersion$statistic, NA, cat_test$homogeneity$`F value`[1]),
    p = c(uniformity$p.value, dispersion$p.value, outliers$p.value, cat_test$homogeneity$`Pr(>F)`[1]),
    row.names = NULL
  )
}

# file-name friendly species label
sp_file <- function(sp) gsub("\\. ", "_", sp)


#Beta regressions
##########################################################################
##########################################################################

lrt_list <- list()
emm_list <- list()
pairs_list <- list()
dharma_list <- list()
model_list <- list()

for (metric in c("AUC", "TSS")) {
  d <- perf %>% filter(metric == !!metric)

  # 1. Overall effect of data source, accounting for species identity
  ##########################################################################
  message("Fitting overall model: ", metric)
  m_species  <- fit_beta(value ~ species, d)
  m_additive <- fit_beta(value ~ source + species, d)
  name_mod <- paste0(metric, "_overall")
  model_list[[name_mod]] <- m_additive
  lrt_list[[name_mod]] <- lrt_row(m_species, m_additive, metric, "1 overall source effect (source + species vs. species)")
  dharma_list[[name_mod]] <- dharma_check(m_additive, d, name_mod)

  # 2. Source x species interaction
  ##########################################################################
  message("Fitting interaction model: ", metric)
  m_interaction <- fit_beta(value ~ source * species, d)
  name_mod <- paste0(metric, "_interaction")
  model_list[[name_mod]] <- m_interaction
  lrt_list[[name_mod]] <- lrt_row(m_additive, m_interaction, metric, "2 source x species interaction (source * species vs. source + species)")
  dharma_list[[name_mod]] <- dharma_check(m_interaction, d, name_mod)

  # 3. Comparison of data sources within each species
  ##########################################################################
  for (sp in levels(perf$species)) {
    message("Fitting species model: ", metric, " - ", sp)
    d_sp <- d %>% filter(species == sp)
    m_null <- fit_beta(value ~ 1, d_sp)
    m_full <- fit_beta(value ~ source, d_sp)
    name_mod <- paste0(metric, "_", sp_file(sp))
    model_list[[name_mod]] <- m_full
    lrt_list[[name_mod]] <- lrt_row(m_null, m_full, metric, "3 source effect within species (source vs. intercept only)", sp)
    dharma_list[[name_mod]] <- dharma_check(m_full, d_sp, name_mod)

    # estimated marginal means (response scale) with 95% CI
    emm <- emmeans(m_full, ~ source, type = "response")
    emm_list[[name_mod]] <- data.frame(metric = metric, species = sp, as.data.frame(confint(emm)))

    # pairwise differences on the response scale, Tukey-adjusted CIs and p-values
    pw <- contrast(regrid(emm), method = "pairwise", adjust = "tukey")
    pairs_list[[name_mod]] <- data.frame(metric = metric, species = sp,
                                         merge(as.data.frame(confint(pw)),
                                               as.data.frame(pw)[, c("contrast", "z.ratio", "p.value")],
                                               by = "contrast"))
  }
}

lrt_table <- bind_rows(lrt_list) %>% arrange(test, metric)
emm_table <- bind_rows(emm_list)
pairs_table <- bind_rows(pairs_list)
dharma_table <- bind_rows(dharma_list)

write.csv(lrt_table, paste0(out_dir, "beta_LRT_tests.csv"), row.names = FALSE)
write.csv(emm_table, paste0(out_dir, "beta_emmeans_within_species.csv"), row.names = FALSE)
write.csv(pairs_table, paste0(out_dir, "beta_pairwise_contrasts_within_species.csv"), row.names = FALSE)
write.csv(dharma_table, paste0(out_dir, "DHARMa_tests.csv"), row.names = FALSE)

print(lrt_table)
print(pairs_table)
print(dharma_table)


#Full model output for documentation
##########################################################################
##########################################################################

sink(paste0(out_dir, "beta_model_summaries.txt"))
for (name_mod in names(model_list)) {
  cat("==========", name_mod, "==========\n")
  print(summary(model_list[[name_mod]]))
  cat("\nLikelihood-ratio test:\n")
  print(lrt_list[[name_mod]])
  if (!is.null(pairs_list[[name_mod]])) {
    cat("\nEstimated marginal means (response scale):\n")
    print(emm_list[[name_mod]])
    cat("\nPairwise contrasts (response scale, Tukey-adjusted):\n")
    print(pairs_list[[name_mod]])
  }
  cat("\nDHARMa residual tests:\n")
  print(dharma_list[[name_mod]])
  cat("\n\n")
}
print(sessionInfo())
sink()


#Supplementary table
##########################################################################
##########################################################################

# data source labels as used in the manuscript
relabel_source <- function(x) {
  x <- gsub("\\bGBIF\\b", "CS", x, perl = TRUE)
  x <- gsub("\\bHLNUG\\b", "StAg", x, perl = TRUE)
  gsub(" - ", " \u2013 ", x)
}
fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))

test_labels <- c(
  "1 overall source effect (source + species vs. species)" = "Overall source effect",
  "2 source x species interaction (source * species vs. source + species)" = "Source \u00d7 species interaction",
  "3 source effect within species (source vs. intercept only)" = "Source effect within species"
)

supp_lrt <- lrt_table %>%
  transmute(
    Test    = test_labels[test],
    Metric  = metric,
    Species = ifelse(species == "all species", "All species", species),
    n       = n_obs,
    chisq   = sprintf("%.2f", chisq),
    df      = df,
    sig     = p < 0.05,
    p       = fmt_p(p)
  )
names(supp_lrt)[names(supp_lrt) == "chisq"] <- "χ²"


supp_pairs <- pairs_table %>%
  arrange(species, metric) %>%
  transmute(
    Species    = species,
    Metric     = metric,
    Contrast   = relabel_source(as.character(contrast)),
    Difference = sprintf("%.3f", estimate),
    `95% CI`   = sprintf("[%.3f, %.3f]", asymp.LCL, asymp.UCL),
    z          = sprintf("%.2f", z.ratio),
    p          = fmt_p(p.value),
    sig        = p.value < 0.05
  )

write.csv(select(supp_lrt, -sig), paste0(out_dir, "Supplementary_table_LRT_tests.csv"), row.names = FALSE)
write.csv(select(supp_pairs, -sig), paste0(out_dir, "Supplementary_table_pairwise_contrasts.csv"), row.names = FALSE)

format_table <- function(tab) {
  ft <- flextable(select(tab, -sig)) %>%
    italic(i = ~ Species != "All species", j = "Species") %>%
    bold(i = tab$sig, j = "p") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "all") %>%
    align(align = "center", part = "all") %>%
    align(j = 1, align = "left", part = "all") %>%
    merge_v(j = 1) %>%
    valign(valign = "top", part = "body") %>%
    hline(border = fp_border(width = 0.5), part = "body") %>%
    autofit()
  fix_border_issues(ft)
}

caption_lrt <- paste(
  "Appendix 4.2: Likelihood-ratio tests of differences in test performance (AUC, TSS) between data sources",
  "(COM, CS, StAg) of the local models across all model runs, based on beta regressions (logit link).",
  "Overall source effect: source + species vs. species only; source \u00d7 species interaction:",
  "source \u00d7 species vs. source + species; source effect within species: source vs. intercept-only model",
  "fitted separately per species (10 runs per data source and species). Significant p-values (< 0.05) in bold."
)
caption_pairs <- paste(
  "Appendix 4.3: Pairwise differences in mean test performance between data sources within each species,",
  "estimated from the species-specific beta regressions (estimated marginal means on the response scale).",
  "Differences are given as first minus second data source; 95% confidence intervals and p-values are",
  "Tukey-adjusted. Significant p-values (< 0.05) in bold."
)
abbreviations <- "COM = combined data; CS = citizen science data (GBIF); StAg = state agency data (HLNUG)."

read_docx() %>%
  body_add_par(caption_lrt) %>%
  body_add_flextable(format_table(supp_lrt)) %>%
  body_add_par(abbreviations) %>%
  body_add_par("") %>%
  body_add_par(caption_pairs) %>%
  body_add_flextable(format_table(supp_pairs)) %>%
  body_add_par(abbreviations) %>%
  print(target = paste0(out_dir, "Supplementary_table_model_comparison.docx"))

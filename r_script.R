# Analysis accompanying:
# Many a many: An experimental investigation on the Polish vague quantifiers
# duzo and wiele
#
############################################################
# MANYDOT 2026
# Polish quantifiers DUZO / WIELE in isolation and together
############################################################

############################
# 0. SETUP
############################

# Run from the directory containing the three input CSV files.
# Required packages (install once if needed):
# install.packages(c("dplyr", "tidyr", "lme4", "ggplot2", "emmeans", "effects"))
# All outputs are written to the current working directory.
input_files <- paste0("results_link_", 1:3, ".csv")
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop("Missing input files: ", paste(missing_files, collapse = ", "))
}

library(dplyr)
library(tidyr)
library(lme4)
library(ggplot2)
library(emmeans)
library(effects)

############################
# 1. LOAD AND COMBINE THE THREE RAW DATASETS
############################

column_names <- c(
  "ResultsReceptionTime",
  "IPHash",
  "Controller",
  "ItemNumber",
  "ElementNumber",
  "Label",
  "LatinSquareGroup",
  "FieldName",
  "FieldValue"
)

duzo <- read.csv(
  "results_link_1.csv",
  header = FALSE,
  comment.char = "#",
  col.names = column_names,
  stringsAsFactors = FALSE
)

wiele <- read.csv(
  "results_link_2.csv",
  header = FALSE,
  comment.char = "#",
  col.names = column_names,
  stringsAsFactors = FALSE
)

both <- read.csv(
  "results_link_3.csv",
  header = FALSE,
  comment.char = "#",
  col.names = column_names,
  stringsAsFactors = FALSE
)

duzo$Dataset <- "DUZO_only"
wiele$Dataset <- "WIELE_only"
both$Dataset <- "BOTH"

raw_data <- rbind(duzo, wiele, both)

############################
# 2. PARTICIPANT INFORMATION AND DESCRIPTIVE STATISTICS
############################

# One record per participant within each dataset.
participants_base <- raw_data %>%
  distinct(Dataset, IPHash) %>%
  mutate(
    Context = ifelse(Dataset == "BOTH", "Together", "Isolation")
  )

# Extract questionnaire data from Form records, if available.
form_raw <- raw_data %>%
  filter(Controller == "Form") %>%
  dplyr::select(Dataset, IPHash, FieldName, FieldValue) %>%
  distinct()

participant_info <- participants_base

if (nrow(form_raw) > 0) {
  form_wide <- form_raw %>%
    tidyr::pivot_wider(names_from = FieldName, values_from = FieldValue)
  
  participant_info <- participants_base %>%
    left_join(form_wide, by = c("Dataset", "IPHash"))
}

# Export the complete participant information table.
write.table(
  participant_info,
  file = "participant_info_raw.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# Participant counts by dataset.
n_participants_dataset <- participant_info %>%
  group_by(Dataset) %>%
  summarise(
    N_participants = n(),
    .groups = "drop"
  )

print(n_participants_dataset)

write.table(
  n_participants_dataset,
  file = "participants_by_dataset.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# Participant counts by context.
n_participants_context <- participant_info %>%
  group_by(Context) %>%
  summarise(
    N_participants = n(),
    .groups = "drop"
  )

print(n_participants_context)

write.table(
  n_participants_context,
  file = "participants_by_context.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# Age descriptives, if age was recorded.
if ("age" %in% names(participant_info)) {
  participant_info$age <- suppressWarnings(as.numeric(participant_info$age))
  
  age_summary_overall <- participant_info %>%
    summarise(
      Mean_age = mean(age, na.rm = TRUE),
      SD_age = sd(age, na.rm = TRUE),
      Min_age = min(age, na.rm = TRUE),
      Max_age = max(age, na.rm = TRUE),
      N_non_missing = sum(!is.na(age))
    )
  
  print(age_summary_overall)
  
  write.table(
    age_summary_overall,
    file = "age_summary_overall.txt",
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  
  age_summary_by_dataset <- participant_info %>%
    group_by(Dataset) %>%
    summarise(
      Mean_age = mean(age, na.rm = TRUE),
      SD_age = sd(age, na.rm = TRUE),
      Min_age = min(age, na.rm = TRUE),
      Max_age = max(age, na.rm = TRUE),
      N_non_missing = sum(!is.na(age)),
      .groups = "drop"
    )
  
  print(age_summary_by_dataset)
  
  write.table(
    age_summary_by_dataset,
    file = "age_summary_by_dataset.txt",
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}

# Helper function to export frequency tables.
save_freq_table <- function(df, varname, filename) {
  if (varname %in% names(df)) {
    tab <- as.data.frame(table(df[[varname]], useNA = "ifany"))
    names(tab) <- c(varname, "Freq")
    print(tab)
    write.table(
      tab,
      file = filename,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  }
}

save_freq_table(participant_info, "region", "freq_region.txt")
save_freq_table(participant_info, "nativelg", "freq_nativelg.txt")
save_freq_table(participant_info, "vision", "freq_vision.txt")
save_freq_table(participant_info, "consent", "freq_consent.txt")
save_freq_table(participant_info, "gender", "freq_gender.txt")
save_freq_table(participant_info, "sex", "freq_sex.txt")

############################
# 3. RETAIN EXPERIMENTAL RESPONSES
############################

response_data <- subset(raw_data, Controller == "MultipleRatings")

############################
# 4. CREATE ANALYSIS VARIABLES
############################

response_data$Answer <- as.numeric(response_data$FieldValue)
response_data$Condition <- response_data$FieldName
response_data$N_Dots <- as.numeric(response_data$ItemNumber) - 1
response_data$Participant <- as.factor(response_data$IPHash)
response_data$Context <- ifelse(response_data$Dataset == "BOTH", "Together", "Isolation")

############################
# 5. RETAIN THE TWO TARGET QUANTIFIERS
############################

target <- subset(response_data, Condition %in% c("ManyDuzo", "ManyWiele"))

############################
# 6. RELABEL CONDITIONS
############################

target$Condition[target$Condition == "ManyDuzo"] <- "DUZO"
target$Condition[target$Condition == "ManyWiele"] <- "WIELE"

############################
# 7. SET VARIABLE TYPES AND FACTOR REFERENCE LEVELS
############################

target$Condition <- factor(target$Condition, levels = c("DUZO", "WIELE"))
target$Context <- factor(target$Context, levels = c("Isolation", "Together"))
target$Participant <- as.factor(target$Participant)
target$N_Dots <- as.numeric(target$N_Dots)
target$Answer <- as.numeric(target$Answer)

############################
# 8. INSPECT AND EXPORT THE TARGET DATASET
############################

str(target)
head(target)
table(target$Condition, target$Context)

write.table(
  target,
  file = "target_DUZO_WIELE_isolation_vs_together.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

############################
# 9. RATING DESCRIPTIVES
############################

# Mean ratings by condition and context.
mean_ratings <- aggregate(Answer ~ Condition + Context, target, mean)
print(mean_ratings)

write.table(
  mean_ratings,
  file = "means_DUZO_WIELE_by_context.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# Mean, standard deviation, and number of observations.
rating_summary <- target %>%
  group_by(Condition, Context) %>%
  summarise(
    Mean = mean(Answer, na.rm = TRUE),
    SD = sd(Answer, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

print(rating_summary)

write.table(
  rating_summary,
  file = "summary_DUZO_WIELE_by_context.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# Descriptive statistics for each number of dots.
ratings_by_numerosity <- target %>%
  group_by(Condition, Context, N_Dots) %>%
  summarise(
    Mean = mean(Answer, na.rm = TRUE),
    SD = sd(Answer, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

print(ratings_by_numerosity)

write.table(
  ratings_by_numerosity,
  file = "DUZO_WIELE_point_by_point.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

############################
# 10. PLOT MEAN RATINGS
############################

mean_plot_data <- target %>%
  group_by(N_Dots, Condition, Context) %>%
  summarise(
    MeanAnswer = mean(Answer, na.rm = TRUE),
    .groups = "drop"
  )

p1 <- ggplot(mean_plot_data, aes(x = N_Dots, y = MeanAnswer, linetype = Condition)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Context) +
  theme_bw() +
  labs(
    title = "DUZO and WIELE: isolation vs. together",
    x = "Number of dots",
    y = "Mean rating"
  )

print(p1)

ggsave(
  "plot_DUZO_WIELE_isolation_vs_together.png",
  plot = p1,
  width = 8,
  height = 5,
  dpi = 300
)

############################
# 11. BOXPLOTS FOR NUMEROSITIES ABOVE 15
############################

target_16 <- subset(target, N_Dots > 15)

target_16$Alternative <- ifelse(target_16$Context == "Together", "Yes", "No")
target_16$Alternative <- factor(target_16$Alternative, levels = c("Yes", "No"))
target_16$N_Dots <- factor(target_16$N_Dots, levels = 16:29)
target_16$Condition <- factor(target_16$Condition, levels = c("DUZO", "WIELE"))

p2 <- ggplot(target_16, aes(x = Answer, y = Alternative, fill = Condition)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.6,
    outlier.size = 1.8
  ) +
  facet_wrap(~ N_Dots, ncol = 1) +
  scale_x_continuous(breaks = 1:5, limits = c(0.8, 5.1)) +
  labs(
    title = "Scores for Duzo and Wiele",
    subtitle = "For number of dots > 15",
    x = "Likert Scale",
    y = "Alternative",
    fill = "Condition"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(p2)

ggsave(
  "Figure_DUZO_WIELE_above15_boxplot.png",
  plot = p2,
  width = 8,
  height = 12,
  dpi = 300
)

############################
# 12. CONDITION-BY-CONTEXT INTERACTION MODEL
############################

model <- lmer(
  Answer ~ Condition * Context + N_Dots + (1 | Participant),
  data = target
)

print(summary(model))

capture.output(
  summary(model),
  file = "model_DUZO_WIELE_isolation_vs_together.txt"
)

############################
# 13. SEPARATE MODELS FOR EACH QUANTIFIER
############################

duzo_only <- subset(target, Condition == "DUZO")
wiele_only <- subset(target, Condition == "WIELE")

model_duzo <- lmer(
  Answer ~ Context + N_Dots + (1 | Participant),
  data = duzo_only
)

model_wiele <- lmer(
  Answer ~ Context + N_Dots + (1 | Participant),
  data = wiele_only
)

print(summary(model_duzo))
print(summary(model_wiele))

capture.output(
  summary(model_duzo),
  file = "model_DUZO_only_context.txt"
)

capture.output(
  summary(model_wiele),
  file = "model_WIELE_only_context.txt"
)

############################
# 14. MODEL COMPARISONS REPORTED IN THE MANUSCRIPT
############################

m0 <- lmer(
  Answer ~ 1 + (1 | Participant),
  data = target,
  REML = FALSE
)

m1 <- lmer(
  Answer ~ N_Dots + (1 | Participant),
  data = target,
  REML = FALSE
)

m2 <- lmer(
  Answer ~ N_Dots + Condition + (1 | Participant),
  data = target,
  REML = FALSE
)

m3 <- lmer(
  Answer ~ N_Dots + Condition + Context + (1 | Participant),
  data = target,
  REML = FALSE
)

m4 <- lmer(
  Answer ~ N_Dots * Condition * Context + (1 | Participant),
  data = target,
  REML = FALSE
)

model_comparison_all <- anova(m0, m1, m2, m3, m4)
print(model_comparison_all)

capture.output(
  model_comparison_all,
  file = "model_comparison_all.txt"
)

model_comparison_m1_m3 <- anova(m1, m3)
print(model_comparison_m1_m3)

capture.output(
  model_comparison_m1_m3,
  file = "model_comparison_m1_m3.txt"
)

############################
# 15. POST-HOC
############################

# Use the additive model m3 for the comparisons reported in the manuscript.
# Because m3 contains no interaction, conditional contrasts are identical
# across levels of the conditioning factor; they are not independent tests.
posthoc_condition <- emmeans(m3, pairwise ~ Condition)
posthoc_context <- emmeans(m3, pairwise ~ Context)
posthoc_condition_by_context <- emmeans(m3, pairwise ~ Condition | Context)
posthoc_context_by_condition <- emmeans(m3, pairwise ~ Context | Condition)

print(posthoc_condition)
print(posthoc_context)
print(posthoc_condition_by_context)
print(posthoc_context_by_condition)

capture.output(
  posthoc_condition,
  file = "posthoc_condition.txt"
)

capture.output(
  posthoc_context,
  file = "posthoc_context.txt"
)

capture.output(
  posthoc_condition_by_context,
  file = "posthoc_condition_by_context.txt"
)

capture.output(
  posthoc_context_by_condition,
  file = "posthoc_context_by_condition.txt"
)

############################
# 16. ALL EFFECTS PLOTS
############################

eff_m3 <- allEffects(m3)

png("all_effects_m3.png", width = 1200, height = 900, res = 150)
plot(eff_m3)
dev.off()

png("effect_condition_context.png", width = 900, height = 700, res = 150)
plot(effect("Condition:Context", m4))
dev.off()

png("effect_ndots_condition.png", width = 900, height = 700, res = 150)
plot(effect("N_Dots:Condition", m4))
dev.off()

############################
# 17. END OF INITIAL ANALYSES
############################

cat("\nInitial analyses completed.\n")
cat("Outputs saved in the current working directory.\n")


############################################################
# MANYDOT 2026
# ADDITIONAL PUBLICATION FIGURES
# 1) Mean ratings across dots
# 2) Model-based marginal means
############################################################

############################
# 3. FIGURE 1: MEAN RATINGS WITH STANDARD ERROR BARS
############################

desc_plot <- target %>%
  group_by(Context, Condition, N_Dots) %>%
  summarise(
    Mean = mean(Answer, na.rm = TRUE),
    SD = sd(Answer, na.rm = TRUE),
    N = n(),
    SE = SD / sqrt(N),
    .groups = "drop"
  )

p1 <- ggplot(desc_plot, aes(x = N_Dots, y = Mean, color = Condition, group = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2, alpha = 0.6) +
  facet_wrap(~ Context, nrow = 1) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 29)) +
  scale_y_continuous(breaks = 1:5) +
  coord_cartesian(ylim = c(1, 5)) +
  labs(
    title = "Acceptability ratings across numerosities",
    x = "Number of dots",
    y = "Mean acceptability rating",
    color = "Condition"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )

print(p1)

ggsave(
  "Figure1_mean_ratings_across_dots.png",
  plot = p1,
  width = 10,
  height = 5,
  dpi = 300
)

############################
# 4. ADDITIVE MODEL USED FOR ESTIMATED MARGINAL MEANS
############################

# Reuse m3 fitted in the model-comparison section above.

############################
# 5. FIGURE 2: ESTIMATED MARGINAL MEANS
############################

emm_ctx_cond <- emmeans(m3, ~ Condition + Context)
emm_df <- as.data.frame(confint(emm_ctx_cond, level = 0.95))
# emmeans names confidence limits differently for finite and infinite df.
if (all(c("lower.CL", "upper.CL") %in% names(emm_df))) {
  emm_df$CI_lower <- emm_df$lower.CL
  emm_df$CI_upper <- emm_df$upper.CL
} else if (all(c("asymp.LCL", "asymp.UCL") %in% names(emm_df))) {
  emm_df$CI_lower <- emm_df$asymp.LCL
  emm_df$CI_upper <- emm_df$asymp.UCL
} else {
  stop("Confidence interval columns were not found in the emmeans output.")
}

p2 <- ggplot(emm_df, aes(x = Context, y = emmean, color = Condition, group = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.08, linewidth = 0.7) +
  scale_y_continuous() +
  labs(
    title = "Model-based estimated marginal means",
    x = "Context",
    y = "Estimated mean rating",
    color = "Condition"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p2)

ggsave(
  "Figure2_model_estimated_marginal_means.png",
  plot = p2,
  width = 7,
  height = 5,
  dpi = 300
)

############################
# 6. EXPORT MANUSCRIPT TABLES
############################

summary_table <- target %>%
  group_by(Condition, Context) %>%
  summarise(
    Mean = mean(Answer, na.rm = TRUE),
    SD = sd(Answer, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

print(summary_table)

write.table(
  summary_table,
  file = "Table_summary_DUZO_WIELE_context.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

write.table(
  emm_df,
  file = "Table_emmeans_condition_context.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

cat("\nFigures and tables created successfully.\n")
cat("Saved:\n")
cat("- Figure1_mean_ratings_across_dots.png\n")
cat("- Figure2_model_estimated_marginal_means.png\n")
cat("- Table_summary_DUZO_WIELE_context.txt\n")
cat("- Table_emmeans_condition_context.txt\n")


# Record the R environment and package versions used for this run.
capture.output(sessionInfo(), file = "session_info.txt")
# Ordinal mixed models for duzo and wiele ratings
# Run this script with the three results_link CSV files in input_dir.
# The Word document explains the coding and interpretation.

# 1. SETUP
# First fit clmm. Change run_bayes to TRUE to also fit brms.
# Install missing packages once, outside this script:
# install.packages(c("ordinal", "brms", "posterior",
#                    "bayesplot", "ggplot2"))
# brms with backend = "rstan" also needs a working C++ toolchain.
input_dir <- "."
out_dir <- "ordinal_analysis_output"
run_bayes <- FALSE
run_prior_check <- FALSE
run_prior_sensitivity <- FALSE
native_codes <- "SLO"
# Corrected sample: retain the primary native-language code only.
# The user identified V2_022 (code "SLO-other") as non-Polish.
# Its 29 ratings are excluded; code "other" remains excluded as before.
# Context coding follows the original analysis and its cell counts.
contexts <- c("Isolation", "Isolation", "Together")
# These IDs are submissions, not verified unique people.
# V2_030 and V2_031 share an IP hash and identical demographics.
# Keep both to match the original data until their status is resolved.
# For example, if V2_031 is a confirmed repeat, use c("V2_031").
exclude_session_ids <- character(0)
set.seed(12345)
if (!requireNamespace("ordinal", quietly = TRUE)) {
  stop("Install ordinal first: install.packages('ordinal')")
}
# Attach ordinal so anova() can resolve its model-comparison method.
library(ordinal)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(x, name) {
  write.csv(x, file.path(out_dir, name), row.names = FALSE)
}

# 2. IMPORT THE THREE RAW EXPORTS
# The files have nine columns, no header, and # comment blocks.
# Use reception time plus IP hash to distinguish submissions.
read_version <- function(path, v, context_label) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  raw <- read.csv(
    path, header = FALSE, comment.char = "#",
    col.names = c("received", "ip_hash", "controller",
                  "raw_item", "element", "label", "group",
                  "field", "value"),
    colClasses = "character", stringsAsFactors = FALSE,
    na.strings = c("", "NULL"), strip.white = TRUE,
    fileEncoding = "UTF-8", check.names = FALSE
  )
  if (anyNA(raw[c("received", "ip_hash", "field", "value")])) {
    stop("Missing submission identifiers or values in ", path)
  }
  raw$session_key <- paste(raw$received, raw$ip_hash, sep = "_")
  keys <- unique(raw$session_key)
  ids <- sprintf("V%d_%03d", v, seq_along(keys))
  get_form_value <- function(field_name) {
    z <- raw[raw$controller == "Form" & raw$field == field_name, ]
    if (anyDuplicated(z$session_key) ||
        !setequal(z$session_key, keys)) {
      stop("Missing or repeated form field: ", field_name)
    }
    z$value[match(keys, z$session_key)]
  }
  meta <- data.frame(
    participant_id = ids, version = paste("Version", v),
    native_code = get_form_value("nativelg"),
    consent = get_form_value("consent"),
    received = raw$received[match(keys, raw$session_key)],
    ip_hash = raw$ip_hash[match(keys, raw$session_key)]
  )
  meta$eligible <- meta$native_code %in% native_codes &
    meta$consent == "yes"
  repeat_hash <- duplicated(meta$ip_hash) |
    duplicated(meta$ip_hash, fromLast = TRUE)
  meta$possible_repeat <- repeat_hash
  meta$explicitly_excluded <- ids %in% exclude_session_ids
  meta$included <- meta$eligible & !meta$explicitly_excluded
  # Export an audit without the IP hash.
  save_csv(meta[setdiff(names(meta), "ip_hash")],
           paste0("session_audit_link_", v, ".csv"))
  if (any(repeat_hash)) {
    warning("Possible repeated submissions in link ", v, ": ",
            paste(ids[repeat_hash], collapse = ", "))
  }
  z <- raw[raw$controller == "MultipleRatings" &
             raw$label == "test" &
             raw$field %in% c("ManyDuzo", "ManyWiele"), ]
  m <- match(z$session_key, keys)
  z <- z[meta$included[m], ]
  m <- match(z$session_key, keys)
  data.frame(
    participant_id = ids[m], version = paste("Version", v),
    context = context_label,
    expression = ifelse(z$field == "ManyDuzo", "duzo", "wiele"),
    item_number = z$raw_item, rating = z$value
  )
}
dat <- do.call(rbind, lapply(seq_len(3), function(v) {
  read_version(file.path(input_dir,
                         paste0("results_link_", v, ".csv")),
               v, contexts[v])
}))
rownames(dat) <- NULL

# 3. PREPARE VARIABLES AND CHECK THE SAMPLE
required <- c("participant_id", "context", "expression",
              "item_number", "rating")
if (!nrow(dat) || anyNA(dat[required])) {
  stop("The analysis data are empty or have missing required values.")
}
if (any(!dat$rating %in% as.character(1:5))) {
  stop("All ratings must be integer categories 1 through 5.")
}
if (anyDuplicated(dat[c("participant_id", "item_number",
                        "expression")])) {
  stop("Repeated rating for the same submission, item and expression.")
}
dat$participant_id <- factor(dat$participant_id)
# Preserve the original shared-item assumption. Check the stimulus key:
# the same item_number must identify the same stimulus across versions.
dat$item_number <- factor(dat$item_number)
dat$context <- factor(dat$context,
                      levels = c("Isolation", "Together"))
dat$expression <- factor(dat$expression, levels = c("duzo", "wiele"))
# Set treatment coding explicitly so the contrast formulas are stable.
contrasts(dat$context) <- contr.treatment(levels(dat$context), base = 1)
contrasts(dat$expression) <- contr.treatment(levels(dat$expression), base = 1)
dat$rating <- ordered(dat$rating, levels = 1:5)
if (any(table(dat$rating) == 0) ||
    any(table(dat$context, dat$expression) == 0)) {
  stop("An empty rating category or design cell needs review.")
}
str(dat)
print(table(dat$rating))
print(table(dat$context, dat$expression))
cat("Included submissions:", nlevels(dat$participant_id), "\n")
# With the supplied files: 110 submissions, 4292 ratings.
# Rating counts 1 to 5: 1181, 484, 480, 652, 1495.
save_csv(dat, "cleaned_duzo_wiele_ratings.csv")
counts <- as.data.frame(table(dat$context, dat$expression, dat$rating))
names(counts) <- c("context", "expression", "rating", "n")
save_csv(counts, "observed_rating_counts.csv")
# Context is tied to version. Its effect may include cohort differences.
# Adding version to this formula cannot separate it from context.

# 4. FIT THE FREQUENTIST ORDINAL MIXED MODEL FIRST
# NEW: clmm fits the ordinal model without Bayesian estimation.
# Both models use a cumulative logit link and flexible thresholds.
# nAGQ = 1 uses Laplace approximation, needed for two random terms.
model_formula <- rating ~ context * expression +
  (1 | participant_id) + (1 | item_number)
fit_clmm <- ordinal::clmm(
  model_formula, data = dat, link = "logit",
  threshold = "flexible", Hess = TRUE, nAGQ = 1
)
print(summary(fit_clmm))
# Inspect convergence warnings, max.grad and cond.H in this summary.
saveRDS(fit_clmm, file.path(out_dir, "fit_clmm.rds"))
# Test the interaction against the nested additive ordinal model.
fit_clmm_additive <- ordinal::clmm(
  rating ~ context + expression +
    (1 | participant_id) + (1 | item_number),
  data = dat, link = "logit", threshold = "flexible",
  Hess = TRUE, nAGQ = 1
)
interaction_test <- anova(fit_clmm_additive, fit_clmm)
print(interaction_test)
capture.output(summary(fit_clmm), interaction_test,
               file = file.path(out_dir, "clmm_summary.txt"))

# 5. PRESERVE THE FIVE ORIGINAL CONTRASTS
# Positive effects indicate a shift toward higher ratings.
coef_names <- c("contextTogether", "expressionwiele",
                "contextTogether:expressionwiele")
if (!all(coef_names %in% names(fit_clmm$beta))) {
  stop("Unexpected coefficient names; check factor coding.")
}
L <- rbind(
  "wiele - duzo: Isolation" = c(0, 1, 0),
  "wiele - duzo: Together" = c(0, 1, 1),
  "Together - Isolation: duzo" = c(1, 0, 0),
  "Together - Isolation: wiele" = c(1, 0, 1),
  "Context x Expression interaction" = c(0, 0, 1)
)
colnames(L) <- coef_names
beta <- fit_clmm$beta[coef_names]
V <- vcov(fit_clmm)[coef_names, coef_names, drop = FALSE]
if (any(!is.finite(V)) || any(eigen(V, symmetric = TRUE)$values <= 0)) {
  stop("Invalid coefficient covariance; resolve clmm convergence.")
}
estimate <- drop(L %*% beta)
se <- sqrt(diag(L %*% V %*% t(L)))
lo <- estimate - qnorm(0.975) * se
hi <- estimate + qnorm(0.975) * se
p <- 2 * pnorm(abs(estimate / se), lower.tail = FALSE)
clmm_contrasts <- data.frame(
  Contrast = rownames(L), Log_odds = estimate, SE = se,
  Lower_95_CI = lo, Upper_95_CI = hi,
  OR = exp(estimate), OR_lower_95 = exp(lo), OR_upper_95 = exp(hi),
  P_unadjusted = p, P_Holm = p.adjust(p, method = "holm"),
  row.names = NULL
)
print(clmm_contrasts)
save_csv(clmm_contrasts, "clmm_contrasts.csv")
# For the interaction, exp(beta) is a ratio of odds ratios.
# Wald CIs are approximate; P_Holm adjusts the five contrast tests.

# 6. OPTIONAL BAYESIAN ORDINAL MODEL IN BRMS
# REVISED: replace MCMCglmm with brms to use pp_check directly.
# The likelihood is ordinal, not Gamma. Gamma is for continuous
# positive responses and does not match bounded 1-5 categories.
# New logit coefficients differ in scale from the old MCMCglmm fit.
if (!run_bayes) {
  message("clmm completed. Set run_bayes <- TRUE for the brms stage.")
} else {
  pkgs <- c("brms", "posterior", "bayesplot", "ggplot2", "rstan")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace,
                               logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop("Install missing packages: ",
         paste(missing_pkgs, collapse = ", "))
  }
  # Explicit regularizing priors on the latent logit scale.
  # normal(0,1) gives most single-coefficient ORs about 0.14 to 7.1.
  # Intercept means ordered cut-points for centered predictors.
  # exponential(1) regularizes nonnegative random-effect SDs.
  priors <- c(
    brms::set_prior("normal(0, 1)", class = "b"),
    brms::set_prior("normal(0, 2.5)", class = "Intercept"),
    brms::set_prior("exponential(1)", class = "sd")
  )
  ncores <- parallel::detectCores()
  if (is.na(ncores)) ncores <- 1L
  fit_bayes <- brms::brm(
    formula = model_formula, data = dat,
    family = brms::cumulative(link = "logit", threshold = "flexible"),
    prior = priors, backend = "rstan",
    chains = 4, cores = min(4L, ncores),
    iter = 4000, warmup = 2000, seed = 12345,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    refresh = 500
  )
  print(summary(fit_bayes))
  print(brms::prior_summary(fit_bayes))
  saveRDS(fit_bayes, file.path(out_dir, "fit_brms.rds"))
  
  # 7. CHECK SAMPLING BEFORE INTERPRETING BAYESIAN RESULTS
  diagnostics <- posterior::summarise_draws(
    posterior::as_draws_array(fit_bayes)
  )
  save_csv(as.data.frame(diagnostics), "brms_diagnostics.csv")
  bad <- with(diagnostics, !is.finite(rhat) | rhat > 1.01 |
                !is.finite(ess_bulk) | ess_bulk < 400 |
                !is.finite(ess_tail) | ess_tail < 400)
  print(diagnostics[bad, ])
  sampler <- rstan::get_sampler_params(fit_bayes$fit,
                                       inc_warmup = FALSE)
  divergences <- sum(vapply(sampler, function(s) {
    sum(s[, "divergent__"])
  }, numeric(1)))
  depth_hits <- sum(vapply(sampler, function(s) {
    sum(s[, "treedepth__"] >= 15)
  }, numeric(1)))
  cat("Divergences:", divergences, "; depth hits:", depth_hits, "\n")
  if (any(bad) || divergences > 0 || depth_hits > 0) {
    warning("Resolve sampling diagnostics before interpreting results.")
  }
  grDevices::pdf(file.path(out_dir, "brms_traceplots.pdf"),
                 width = 11, height = 8)
  plot(fit_bayes, variable = "^b_", regex = TRUE, ask = FALSE)
  grDevices::dev.off()
  
  # 8. POSTERIOR PREDICTIVE CHECKS AGAINST THE OBSERVED RATINGS
  # NEW: bar plots respect discrete categories; use no density smoothing.
  # These are in-sample checks including participant and item effects.
  ppc <- list(
    overall = brms::pp_check(fit_bayes, type = "bars",
                             ndraws = 200, prob = 0.95, freq = FALSE),
    context = brms::pp_check(fit_bayes, type = "bars_grouped",
                             group = "context", ndraws = 200,
                             prob = 0.95, freq = FALSE),
    expression = brms::pp_check(fit_bayes, type = "bars_grouped",
                                group = "expression", ndraws = 200,
                                prob = 0.95, freq = FALSE)
  )
  # Check all four context-by-expression cells, not just marginal fit.
  yrep <- brms::posterior_predict(fit_bayes, ndraws = 200)
  cells <- interaction(dat$context, dat$expression, sep = " / ")
  ppc$cells <- bayesplot::ppc_bars_grouped(
    y = as.integer(dat$rating), yrep = yrep, group = cells,
    prob = 0.95, freq = FALSE
  )
  for (nm in names(ppc)) {
    print(ppc[[nm]])
    ggplot2::ggsave(file.path(out_dir, paste0("ppc_", nm, ".png")),
                    plot = ppc[[nm]], width = 9, height = 6, dpi = 200)
  }
  # Look for systematic misses, especially categories 1 and 5 and
  # individual cells. PPCs do not prove proportional odds or held-out fit.
  
  # 9. BAYESIAN VERSIONS OF THE SAME FIVE CONTRASTS
  beta_draws <- brms::fixef(fit_bayes, summary = FALSE)
  if (!all(coef_names %in% colnames(beta_draws))) {
    stop("Unexpected brms coefficient names; check factor coding.")
  }
  cd <- beta_draws[, coef_names, drop = FALSE] %*% t(L)
  summarize_contrast <- function(x) {
    q <- unname(quantile(x, c(0.025, 0.975)))
    c(Mean_log_odds = mean(x), Median_log_odds = median(x),
      Lower_95_CrI = q[1], Upper_95_CrI = q[2],
      Pr_greater_zero = mean(x > 0), Pr_less_zero = mean(x < 0),
      Median_OR = exp(median(x)),
      OR_lower_95 = exp(q[1]), OR_upper_95 = exp(q[2]))
  }
  tab <- t(apply(cd, 2, summarize_contrast))
  bayes_contrasts <- data.frame(Contrast = rownames(tab), tab,
                                row.names = NULL, check.names = FALSE)
  print(bayes_contrasts)
  save_csv(bayes_contrasts, "brms_contrasts.csv")
  # CrIs here are equal-tailed, not HPD intervals. Directional posterior
  # probabilities are not p-values or probabilities of the null model.
  
  # 10. OPTIONAL PRIOR CHECK AND PRIOR SENSITIVITY
  # Run these before final reporting to assess the chosen prior scales.
  if (run_prior_check) {
    fit_prior <- update(fit_bayes, sample_prior = "only",
                        iter = 2000, warmup = 1000, seed = 12346)
    prior_plot <- brms::pp_check(fit_prior, type = "bars",
                                 ndraws = 200, prob = 0.95, freq = FALSE)
    print(prior_plot)
    ggplot2::ggsave(file.path(out_dir, "prior_predictive.png"),
                    prior_plot, width = 9, height = 6, dpi = 200)
  }
  if (run_prior_sensitivity) {
    wider_priors <- c(
      brms::set_prior("normal(0, 2)", class = "b"),
      brms::set_prior("normal(0, 2.5)", class = "Intercept"),
      brms::set_prior("exponential(0.5)", class = "sd")
    )
    fit_wider <- update(fit_bayes, prior = wider_priors, seed = 12347)
    print(summary(fit_wider))
    saveRDS(fit_wider, file.path(out_dir, "fit_brms_wider_priors.rds"))
    wd <- brms::fixef(fit_wider, summary = FALSE)
    wc <- wd[, coef_names, drop = FALSE] %*% t(L)
    wt <- t(apply(wc, 2, summarize_contrast))
    save_csv(data.frame(Contrast = rownames(wt), wt, row.names = NULL),
             "brms_contrasts_wider_priors.csv")
    # Compare contrast tables; repeat sampling and PPC checks for fit_wider.
  }
}
capture.output(sessionInfo(), file = file.path(out_dir, "session_info.txt"))
# Both cumulative-logit models assume proportional odds: each predictor
# has the same slope across thresholds. Assess this before final reporting.
# If fit is poor, review threshold-specific effects or a different ordinal
# family; poor fit alone does not justify a continuous Gamma likelihood.

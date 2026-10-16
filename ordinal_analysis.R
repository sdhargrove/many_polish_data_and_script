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

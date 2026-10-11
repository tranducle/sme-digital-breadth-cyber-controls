suppressPackageStartupMessages(library(survey))

if (as.character(getRversion()) != "4.6.1") stop("R version drift")
if (as.character(packageVersion("survey")) != "4.5") stop("survey version drift")

uk_path <- "results/uk/build/uk_survey_design.rds"
it_path <- "results/italy/build/italy_survey_design.rds"
out_dir <- "results/size_gap"

uk_des <- readRDS(uk_path)
it_des <- readRDS(it_path)

if (nrow(uk_des$variables) != 2859) stop("UK N drift")
if (nrow(it_des$variables) != 16809) stop("Italy N drift")

numeric_gradient <- function(fn, beta_value) {
  g <- numeric(length(beta_value))
  names(g) <- names(beta_value)
  for (j in seq_along(beta_value)) {
    h <- 1e-6 * max(1, abs(beta_value[j]))
    bp <- beta_value
    bm <- beta_value
    bp[j] <- bp[j] + h
    bm[j] <- bm[j] - h
    g[j] <- (fn(bp) - fn(bm)) / (2 * h)
  }
  g
}

delta_contrast <- function(model, data_obj, weights_vec, set_small, set_medium) {
  beta <- coef(model)
  vc <- vcov(model)
  trm <- delete.response(terms(model))

  standardized_mean <- function(beta_value, state_fn) {
    nd <- data_obj
    nd <- state_fn(nd)
    x <- model.matrix(trm, data=nd)
    x <- x[, names(beta), drop=FALSE]
    pr <- plogis(drop(x %*% beta_value))
    sum(weights_vec * pr) / sum(weights_vec)
  }

  contrast_fn <- function(beta_value) {
    standardized_mean(beta_value, set_medium) - standardized_mean(beta_value, set_small)
  }

  est <- contrast_fn(beta)
  grad <- numeric_gradient(contrast_fn, beta)
  se <- sqrt(drop(t(grad) %*% vc %*% grad))
  crit <- qt(0.975, df=model$df.residual)
  c(
    estimate=est,
    se=se,
    lower=est - crit * se,
    upper=est + crit * se
  )
}

make_row <- function(country, spec, stage, model, data_obj, weights_vec, set_small, set_medium) {
  dr <- delta_contrast(model, data_obj, weights_vec, set_small, set_medium)
  data.frame(
    country=country,
    specification=spec,
    stage=stage,
    estimate=unname(dr["estimate"]),
    se=unname(dr["se"]),
    lower=unname(dr["lower"]),
    upper=unname(dr["upper"]),
    model_df=model$df.residual,
    converged=isTRUE(model$converged)
  )
}

run_country <- function(country, design_obj, weight_name, size_mode, base_formula, conditioned_formula, categorical_formula) {
  d <- design_obj$variables
  w <- d[[weight_name]]

  if (size_mode == "uk") {
    set_small <- function(x) { x$size_medium <- 0L; x }
    set_medium <- function(x) { x$size_medium <- 1L; x }
  } else {
    size_levels <- levels(d$size_factor)
    set_small <- function(x) { x$size_factor <- factor("Small", levels=size_levels); x }
    set_medium <- function(x) { x$size_factor <- factor("Medium", levels=size_levels); x }
  }

  m_abs <- svyglm(base_formula, design=design_obj, family=quasibinomial(link="logit"), rescale=FALSE)
  m_cond <- svyglm(conditioned_formula, design=design_obj, family=quasibinomial(link="logit"), rescale=FALSE)

  primary <- rbind(
    make_row(country, "primary_weighted_continuous", "absolute", m_abs, d, w, set_small, set_medium),
    make_row(country, "primary_weighted_continuous", "conditioned", m_cond, d, w, set_small, set_medium)
  )

  d_eq <- d
  d_eq$equal_weight <- 1
  des_eq <- svydesign(ids=~1, data=d_eq, weights=~equal_weight)

  m_abs_eq <- svyglm(base_formula, design=des_eq, family=quasibinomial(link="logit"), rescale=FALSE)
  m_cond_eq <- svyglm(conditioned_formula, design=des_eq, family=quasibinomial(link="logit"), rescale=FALSE)

  equal_weight <- rbind(
    make_row(country, "S1_equal_weight_continuous", "absolute", m_abs_eq, d_eq, d_eq$equal_weight, set_small, set_medium),
    make_row(country, "S1_equal_weight_continuous", "conditioned", m_cond_eq, d_eq, d_eq$equal_weight, set_small, set_medium)
  )

  des_cat <- update(design_obj, footprint_factor=factor(footprint_score, levels=0:4))
  d_cat <- des_cat$variables
  w_cat <- d_cat[[weight_name]]

  m_abs_cat <- svyglm(base_formula, design=des_cat, family=quasibinomial(link="logit"), rescale=FALSE)
  m_cond_cat <- svyglm(categorical_formula, design=des_cat, family=quasibinomial(link="logit"), rescale=FALSE)

  categorical <- rbind(
    make_row(country, "S2_weighted_categorical_footprint", "absolute", m_abs_cat, d_cat, w_cat, set_small, set_medium),
    make_row(country, "S2_weighted_categorical_footprint", "conditioned", m_cond_cat, d_cat, w_cat, set_small, set_medium)
  )

  rows <- rbind(primary, equal_weight, categorical)
  list(
    rows=rows,
    models=list(
      primary_absolute=m_abs,
      primary_conditioned=m_cond,
      equal_absolute=m_abs_eq,
      equal_conditioned=m_cond_eq,
      categorical_absolute=m_abs_cat,
      categorical_conditioned=m_cond_cat
    )
  )
}

uk <- run_country(
  country="UK",
  design_obj=uk_des,
  weight_name="w_norm",
  size_mode="uk",
  base_formula=response_breadth ~ size_medium + year_factor + region_factor,
  conditioned_formula=response_breadth ~ size_medium + footprint_breadth + year_factor + region_factor,
  categorical_formula=response_breadth ~ size_medium + footprint_factor + year_factor + region_factor
)

it <- run_country(
  country="Italy",
  design_obj=it_des,
  weight_name="weight",
  size_mode="italy",
  base_formula=response_breadth ~ size_factor + region_factor + macrosector_factor,
  conditioned_formula=response_breadth ~ size_factor + footprint_breadth + region_factor + macrosector_factor,
  categorical_formula=response_breadth ~ size_factor + footprint_factor + region_factor + macrosector_factor
)

contrasts <- rbind(uk$rows, it$rows)

make_attenuation <- function(df) {
  specs <- unique(df$specification)
  out <- list()
  for (sp in specs) {
    z <- df[df$specification == sp, ]
    abs_est <- z$estimate[z$stage == "absolute"]
    cond_est <- z$estimate[z$stage == "conditioned"]
    if (length(abs_est) != 1 || length(cond_est) != 1) stop("contrast extraction failure")
    attenuation <- if (abs(abs_est) > 0) 1 - abs(cond_est) / abs(abs_est) else NA_real_
    out[[length(out)+1]] <- data.frame(
      country=unique(z$country),
      specification=sp,
      absolute_contrast=abs_est,
      conditioned_contrast=cond_est,
      signed_change=cond_est - abs_est,
      attenuation_fraction=attenuation,
      attenuation_percent=100 * attenuation,
      same_sign=sign(abs_est) == sign(cond_est),
      conditioned_abs_le_absolute_abs=abs(cond_est) <= abs(abs_est)
    )
  }
  do.call(rbind, out)
}

att <- rbind(
  make_attenuation(contrasts[contrasts$country == "UK", ]),
  make_attenuation(contrasts[contrasts$country == "Italy", ])
)

primary_att <- att[att$specification == "primary_weighted_continuous", ]
uk_att <- primary_att$attenuation_fraction[primary_att$country == "UK"]
it_att <- primary_att$attenuation_fraction[primary_att$country == "Italy"]
uk_no_increase <- primary_att$conditioned_abs_le_absolute_abs[primary_att$country == "UK"]
it_no_increase <- primary_att$conditioned_abs_le_absolute_abs[primary_att$country == "Italy"]
same_sign_both <- all(primary_att$same_sign)

if (
  same_sign_both &&
  isTRUE(uk_att >= 0.20) &&
  isTRUE(it_att >= 0.20) &&
  isTRUE(uk_no_increase) &&
  isTRUE(it_no_increase)
) {
  primary_class <- "STRONG"
} else if (
  same_sign_both &&
  isTRUE(uk_att > 0) &&
  isTRUE(it_att > 0) &&
  isTRUE(uk_no_increase) &&
  isTRUE(it_no_increase) &&
  (isTRUE(uk_att >= 0.20) || isTRUE(it_att >= 0.20))
) {
  primary_class <- "MODERATE"
} else {
  primary_class <- "WEAK_OR_ABSENT"
}

sens <- att[att$specification != "primary_weighted_continuous", ]
sensitivity_reversal <- any(!sens$conditioned_abs_le_absolute_abs)
sensitivity_all_attenuate <- all(sens$attenuation_fraction > 0)

checks <- c(
  uk_n_match=nrow(uk_des$variables) == 2859,
  italy_n_match=nrow(it_des$variables) == 16809,
  all_models_converged=all(contrasts$converged),
  all_contrasts_finite=all(is.finite(contrasts$estimate)),
  all_standard_errors_finite_positive=all(is.finite(contrasts$se)) && all(contrasts$se > 0),
  all_confidence_limits_finite=all(is.finite(contrasts$lower)) && all(is.finite(contrasts$upper)),
  attenuation_finite=all(is.finite(att$attenuation_fraction)),
  no_primary_sign_flip=all(primary_att$same_sign),
  sensitivity_no_reversal=!sensitivity_reversal
)

write.csv(contrasts, file=file.path(out_dir, "size_gap_size_contrasts.csv"), row.names=FALSE)
write.csv(att, file=file.path(out_dir, "size_gap_attenuation_summary.csv"), row.names=FALSE)
writeLines(paste(names(checks), checks, sep="="), file.path(out_dir, "size_gap_specification_checks.txt"))
writeLines(c(
  paste("primary_classification", primary_class, sep="="),
  paste("uk_primary_attenuation_percent", 100 * uk_att, sep="="),
  paste("italy_primary_attenuation_percent", 100 * it_att, sep="="),
  paste("sensitivity_all_attenuate", sensitivity_all_attenuate, sep="="),
  paste("sensitivity_reversal", sensitivity_reversal, sep="="),
  paste("r_version", as.character(getRversion()), sep="="),
  paste("survey_version", as.character(packageVersion("survey")), sep="=")
), file.path(out_dir, "size_gap_decision.txt"))

saveRDS(
  list(uk=uk$models, italy=it$models),
  file.path(out_dir, "size_gap_models.rds"),
  version=3
)

if (!all(checks)) quit(status=2)

cat("SIZE_GAP_EXECUTION_PASS\n")
cat("primary_classification=", primary_class, "\n", sep="")
print(primary_att)
cat("sensitivities\n")
print(sens)

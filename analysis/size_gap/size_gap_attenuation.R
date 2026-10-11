suppressPackageStartupMessages(library(survey))

if (as.character(getRversion()) != "4.6.1") stop("R version drift")
if (as.character(packageVersion("survey")) != "4.5") stop("survey version drift")

root <- "."
out_dir <- "results/size_gap"

uk_weighted <- readRDS(file.path(root, "results/uk/build/uk_survey_design.rds"))
it_weighted <- readRDS(file.path(root, "results/italy/build/italy_survey_design.rds"))

if (nrow(uk_weighted$variables) != 2859) stop("UK N drift")
if (nrow(it_weighted$variables) != 16809) stop("Italy N drift")

std_contrast <- function(model, data_obj, weights, size_var, small_value, medium_value) {
  beta <- coef(model)
  vc <- vcov(model)
  terms_no_y <- delete.response(terms(model))
  df_resid <- model$df.residual

  predict_avg <- function(beta_value, size_value) {
    nd <- data_obj
    if (is.factor(nd[[size_var]])) {
      nd[[size_var]] <- factor(size_value, levels=levels(data_obj[[size_var]]))
    } else {
      nd[[size_var]] <- size_value
    }
    X <- model.matrix(terms_no_y, data=nd)
    X <- X[, names(beta), drop=FALSE]
    p <- plogis(drop(X %*% beta_value))
    sum(weights * p) / sum(weights)
  }

  contrast_fn <- function(beta_value) {
    predict_avg(beta_value, medium_value) - predict_avg(beta_value, small_value)
  }

  est <- contrast_fn(beta)
  g <- numeric(length(beta))
  names(g) <- names(beta)
  for (j in seq_along(beta)) {
    h <- 1e-6 * max(1, abs(beta[j]))
    bp <- beta
    bm <- beta
    bp[j] <- bp[j] + h
    bm[j] <- bm[j] - h
    g[j] <- (contrast_fn(bp) - contrast_fn(bm)) / (2 * h)
  }
  se <- sqrt(drop(t(g) %*% vc %*% g))
  crit <- qt(0.975, df=df_resid)
  c(estimate=est, se=se, lower=est-crit*se, upper=est+crit*se)
}

run_pair <- function(country, design, formula_abs, formula_cond, size_var, small_value, medium_value, weight_name, analysis_name) {
  d <- design$variables
  w <- d[[weight_name]]
  ma <- svyglm(formula_abs, design=design, family=quasibinomial(link="logit"), rescale=FALSE)
  mc <- svyglm(formula_cond, design=design, family=quasibinomial(link="logit"), rescale=FALSE)
  if (!isTRUE(ma$converged) || !isTRUE(mc$converged)) stop(paste(country, analysis_name, "convergence failure"))
  a <- std_contrast(ma, d, w, size_var, small_value, medium_value)
  cnd <- std_contrast(mc, d, w, size_var, small_value, medium_value)
  attenuation <- if (abs(a["estimate"]) > 1e-12) 1 - abs(cnd["estimate"]) / abs(a["estimate"]) else NA_real_
  signed_change <- cnd["estimate"] - a["estimate"]
  data.frame(
    country=country,
    analysis=analysis_name,
    absolute_contrast=unname(a["estimate"]),
    absolute_se=unname(a["se"]),
    absolute_lower=unname(a["lower"]),
    absolute_upper=unname(a["upper"]),
    conditioned_contrast=unname(cnd["estimate"]),
    conditioned_se=unname(cnd["se"]),
    conditioned_lower=unname(cnd["lower"]),
    conditioned_upper=unname(cnd["upper"]),
    attenuation=unname(attenuation),
    signed_change=unname(signed_change),
    same_sign=(sign(a["estimate"]) == sign(cnd["estimate"]))
  )
}

make_equal_design <- function(design, weight_name) {
  d <- design$variables
  d$one_weight <- 1
  svydesign(ids=~1, weights=~one_weight, data=d)
}

uk_abs <- response_breadth ~ size_medium + year_factor + region_factor
uk_cond <- response_breadth ~ footprint_breadth + size_medium + year_factor + region_factor
it_abs <- response_breadth ~ size_factor + region_factor + macrosector_factor
it_cond <- response_breadth ~ footprint_breadth + size_factor + region_factor + macrosector_factor

primary <- rbind(
  run_pair("UK", uk_weighted, uk_abs, uk_cond, "size_medium", 0, 1, "w_norm", "primary_weighted_continuous"),
  run_pair("Italy", it_weighted, it_abs, it_cond, "size_factor", "Small", "Medium", "weight", "primary_weighted_continuous")
)

uk_equal <- make_equal_design(uk_weighted, "w_norm")
it_equal <- make_equal_design(it_weighted, "weight")

s1 <- rbind(
  run_pair("UK", uk_equal, uk_abs, uk_cond, "size_medium", 0, 1, "one_weight", "S1_equal_weight"),
  run_pair("Italy", it_equal, it_abs, it_cond, "size_factor", "Small", "Medium", "one_weight", "S1_equal_weight")
)

uk_cat <- update(uk_weighted, footprint_factor=factor(footprint_score, levels=0:4))
it_cat <- update(it_weighted, footprint_factor=factor(footprint_score, levels=0:4))

uk_cat_abs <- response_breadth ~ size_medium + year_factor + region_factor
uk_cat_cond <- response_breadth ~ footprint_factor + size_medium + year_factor + region_factor
it_cat_abs <- response_breadth ~ size_factor + region_factor + macrosector_factor
it_cat_cond <- response_breadth ~ footprint_factor + size_factor + region_factor + macrosector_factor

s2 <- rbind(
  run_pair("UK", uk_cat, uk_cat_abs, uk_cat_cond, "size_medium", 0, 1, "w_norm", "S2_categorical_footprint"),
  run_pair("Italy", it_cat, it_cat_abs, it_cat_cond, "size_factor", "Small", "Medium", "weight", "S2_categorical_footprint")
)

all_results <- rbind(primary, s1, s2)
write.csv(all_results, file.path(out_dir, "size_gap_all_specifications.csv"), row.names=FALSE)

p <- subset(all_results, analysis=="primary_weighted_continuous")
if (nrow(p) != 2) stop("primary row count drift")

strong <- all(p$same_sign) &&
          all(p$attenuation >= 0.20) &&
          all(abs(p$conditioned_contrast) <= abs(p$absolute_contrast))

moderate <- !strong &&
            any(p$attenuation >= 0.20) &&
            all(p$attenuation > 0) &&
            all(abs(p$conditioned_contrast) <= abs(p$absolute_contrast))

sens <- subset(all_results, analysis!="primary_weighted_continuous")
sens_no_reversal <- all(sens$attenuation > 0)

verdict <- if (strong && sens_no_reversal) {
  "STRONG_INFERENTIAL_CONSEQUENCE"
} else if ((strong || moderate) && sens_no_reversal) {
  "MODERATE_INFERENTIAL_CONSEQUENCE"
} else {
  "WEAK_OR_ABSENT_INFERENTIAL_CONSEQUENCE"
}

decision <- data.frame(
  verdict=verdict,
  uk_primary_attenuation=p$attenuation[p$country=="UK"],
  italy_primary_attenuation=p$attenuation[p$country=="Italy"],
  all_primary_same_sign=all(p$same_sign),
  all_primary_nonincrease=all(abs(p$conditioned_contrast) <= abs(p$absolute_contrast)),
  sensitivity_attenuation_no_reversal=sens_no_reversal
)
write.csv(decision, file.path(out_dir, "size_gap_decision.csv"), row.names=FALSE)

checks <- c(
  uk_n=nrow(uk_weighted$variables)==2859,
  italy_n=nrow(it_weighted$variables)==16809,
  finite_results=all(is.finite(unlist(all_results[,c(
    "absolute_contrast","absolute_se","absolute_lower","absolute_upper",
    "conditioned_contrast","conditioned_se","conditioned_lower","conditioned_upper",
    "attenuation","signed_change"
  )]))),
  positive_se=all(all_results$absolute_se > 0) && all(all_results$conditioned_se > 0),
  primary_rows=nrow(p)==2,
  sensitivity_rows=nrow(sens)==4
)
writeLines(paste(names(checks), checks, sep="="), file.path(out_dir, "size_gap_attenuation_checks.txt"))
writeLines(c(
  paste("r_version", as.character(getRversion()), sep="="),
  paste("survey_version", as.character(packageVersion("survey")), sep="="),
  paste("verdict", verdict, sep="=")
), file.path(out_dir, "size_gap_manifest.txt"))

if (!all(checks)) quit(status=2)

cat("SIZE_GAP_EXECUTION_PASS\n")
print(all_results)
print(decision)

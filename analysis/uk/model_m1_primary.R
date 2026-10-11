suppressPackageStartupMessages(library(survey))

if (as.character(getRversion()) != "4.6.1") stop("R version drift")
if (as.character(packageVersion("survey")) != "4.5") stop("survey version drift")

design <- readRDS("results/uk/build/uk_survey_design.rds")
d <- design$variables

if (nrow(d) != 2859) stop("analysis N drift")
if (!all(d$manage1 %in% c(0,1))) stop("manage1 coding drift")

f1 <- response_breadth ~ footprint_breadth * manage1 + year_factor + size_medium + region_factor

m1 <- svyglm(
  f1,
  design=design,
  family=quasibinomial(link="logit"),
  rescale=FALSE
)

if (!isTRUE(m1$converged)) stop("M1 did not converge")

coef_tab <- summary(m1)$coefficients
vc <- vcov(m1)
beta <- coef(m1)
df_resid <- m1$df.residual
w <- d$w_norm
interaction_term <- "footprint_breadth:manage1"

make_x <- function(data_obj) {
  x <- model.matrix(delete.response(terms(f1)), data=data_obj)
  x[, names(beta), drop=FALSE]
}

weighted_prediction <- function(beta_value, data_obj) {
  x <- make_x(data_obj)
  p <- plogis(drop(x %*% beta_value))
  sum(w * p) / sum(w)
}

weighted_ame_state <- function(beta_value, state_value) {
  nd <- d
  nd$manage1 <- state_value
  x <- make_x(nd)
  p <- plogis(drop(x %*% beta_value))
  slope <- unname(beta_value["footprint_breadth"]) +
    state_value * unname(beta_value[interaction_term])
  sum(w * slope * p * (1 - p)) / sum(w)
}

num_grad <- function(fn, beta_value) {
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

delta_result <- function(fn, beta_value, vcov_mat, df_resid) {
  est <- fn(beta_value)
  g <- num_grad(fn, beta_value)
  se <- sqrt(drop(t(g) %*% vcov_mat %*% g))
  crit <- qt(0.975, df=df_resid)
  c(estimate=est, se=se, lower=est-crit*se, upper=est+crit*se)
}

ame0_fn <- function(b) weighted_ame_state(b, 0)
ame1_fn <- function(b) weighted_ame_state(b, 1)
amediff_fn <- function(b) weighted_ame_state(b, 1) - weighted_ame_state(b, 0)

ame0 <- delta_result(ame0_fn, beta, vc, df_resid)
ame1 <- delta_result(ame1_fn, beta, vc, df_resid)
amediff <- delta_result(amediff_fn, beta, vc, df_resid)

ame_out <- rbind(
  data.frame(state=0, estimand="weighted_average_marginal_effect", estimate=unname(ame0["estimate"]), se=unname(ame0["se"]), lower=unname(ame0["lower"]), upper=unname(ame0["upper"])),
  data.frame(state=1, estimand="weighted_average_marginal_effect", estimate=unname(ame1["estimate"]), se=unname(ame1["se"]), lower=unname(ame1["lower"]), upper=unname(ame1["upper"]))
)

diff_out <- data.frame(
  contrast="AME_manage1_1_minus_0",
  estimate=unname(amediff["estimate"]),
  se=unname(amediff["se"]),
  df=df_resid,
  lower=unname(amediff["lower"]),
  upper=unname(amediff["upper"])
)

pred_rows <- list()
for (state in c(0,1)) {
  for (score in 0:4) {
    nd <- d
    nd$manage1 <- state
    nd$footprint_breadth <- score / 4
    fn <- function(b) weighted_prediction(b, nd)
    dr <- delta_result(fn, beta, vc, df_resid)
    pred_rows[[length(pred_rows)+1]] <- data.frame(
      manage1=state,
      footprint_score=score,
      footprint_breadth=score/4,
      prediction=unname(dr["estimate"]),
      se=unname(dr["se"]),
      lower=unname(dr["lower"]),
      upper=unname(dr["upper"])
    )
  }
}
pred_out <- do.call(rbind, pred_rows)

est <- coef_tab[interaction_term, 1]
se <- coef_tab[interaction_term, 2]
pval <- coef_tab[interaction_term, 4]
crit <- qt(0.975, df=df_resid)
primary <- data.frame(
  term=interaction_term,
  estimate=est,
  se=se,
  df=df_resid,
  p_value=pval,
  lower=est-crit*se,
  upper=est+crit*se
)

fit <- fitted(m1)
checks <- c(
  version_match=TRUE,
  n_match=nrow(d)==2859,
  converged=isTRUE(m1$converged),
  finite_coefficients=all(is.finite(coef_tab[,1])),
  finite_positive_se=all(is.finite(coef_tab[,2])) && all(coef_tab[,2] > 0),
  covariance_valid=all(is.finite(diag(vc))) && all(diag(vc) >= 0),
  fitted_unit_interval=all(is.finite(fit)) && all(fit >= 0 & fit <= 1),
  interaction_finite=all(is.finite(unlist(primary[,c("estimate","se","p_value","lower","upper")]))),
  ame_states_finite=all(is.finite(unlist(ame_out[,c("estimate","se","lower","upper")]))),
  ame_difference_finite=all(is.finite(unlist(diff_out[,c("estimate","se","lower","upper")]))),
  predictions_finite=all(is.finite(unlist(pred_out[,c("prediction","se","lower","upper")]))),
  predictions_unit_interval=all(pred_out$prediction >= 0 & pred_out$prediction <= 1)
)

write.csv(primary, file=file.path("results/uk/primary", "M1_interaction_coefficient.csv"), row.names=FALSE)
write.csv(ame_out, file=file.path("results/uk/primary", "M1_state_ames.csv"), row.names=FALSE)
write.csv(diff_out, file=file.path("results/uk/primary", "M1_ame_difference.csv"), row.names=FALSE)
write.csv(pred_out, file=file.path("results/uk/primary", "M1_standardized_predictions.csv"), row.names=FALSE)
write.csv(data.frame(term=rownames(coef_tab), estimate=coef_tab[,1], se=coef_tab[,2], statistic=coef_tab[,3], p_value=coef_tab[,4], row.names=NULL), file=file.path("results/uk/primary", "M1_full_coefficients.csv"), row.names=FALSE)
writeLines(paste(names(checks), checks, sep="="), file.path("results/uk/primary", "M1_integrity_checks.txt"))
writeLines(c(
  paste("n", nrow(d), sep="="),
  paste("model_df_residual", df_resid, sep="="),
  paste("converged", m1$converged, sep="="),
  paste("dispersion", summary(m1)$dispersion, sep="="),
  paste("r_version", as.character(getRversion()), sep="="),
  paste("survey_version", as.character(packageVersion("survey")), sep="=")
), file.path("results/uk/primary", "M1_model_manifest.txt"))

saveRDS(m1, file.path("results/uk/primary", "M1_model.rds"), version=3)

if (!all(checks)) quit(status=2)

cat("M1_EXECUTION_PASS\n")
cat("interaction_estimate=", est, "\n", sep="")
cat("interaction_p=", pval, "\n", sep="")
cat("ame0=", ame_out$estimate[1], "\n", sep="")
cat("ame1=", ame_out$estimate[2], "\n", sep="")
cat("ame_difference=", diff_out$estimate, "\n", sep="")

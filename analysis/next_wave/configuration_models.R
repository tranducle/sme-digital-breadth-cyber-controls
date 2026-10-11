options(stringsAsFactors=FALSE)

input_dir <- "results/next_wave/inputs"
out_dir <- "results/next_wave/configuration"

adj <- read.csv(file.path(input_dir, "configuration_adjacent_pairs.csv"), check.names=FALSE)
gap2 <- read.csv(file.path(input_dir, "configuration_two_wave_gap_pairs.csv"), check.names=FALSE)
elig <- read.csv(file.path(input_dir, "configuration_retention_eligibility.csv"), check.names=FALSE)

if (nrow(adj) != 447) stop("adjacent N drift")
if (length(unique(adj$id)) != 308) stop("adjacent organization drift")
if (sum(adj$I13_next) != 268) stop("primary event drift")
if (nrow(gap2) != 154) stop("gap2 N drift")
if (length(unique(gap2$id)) != 116) stop("gap2 organization drift")
if (sum(gap2$I13_next) != 105) stop("gap2 event drift")

cluster_vcov <- function(model, cluster) {
  X <- model.matrix(model)
  y <- model.response(model.frame(model))
  mu <- fitted(model)
  pw <- model$prior.weights
  if (is.null(pw)) pw <- rep(1, length(y))
  if (length(cluster) != length(y)) stop("cluster length mismatch")

  score_rows <- X * as.numeric(pw * (y - mu))
  score_cluster <- rowsum(score_rows, group=cluster, reorder=FALSE)
  meat <- crossprod(score_cluster)

  info_weights <- as.numeric(pw * mu * (1 - mu))
  info <- crossprod(X, X * info_weights)
  if (qr(info)$rank != ncol(info)) stop("information matrix rank failure")
  bread <- solve(info)

  N <- nrow(X)
  K <- ncol(X)
  G <- length(unique(cluster))
  if (G <= 1 || N <= K) stop("invalid CR1 dimensions")
  cr1 <- (G / (G - 1)) * ((N - 1) / (N - K))
  V <- cr1 * bread %*% meat %*% bread

  if (!all(is.finite(V))) stop("nonfinite cluster covariance")
  list(V=V, G=G, N=N, K=K, cr1=cr1)
}

num_grad <- function(fn, beta) {
  g <- numeric(length(beta))
  names(g) <- names(beta)
  for (j in seq_along(beta)) {
    h <- 1e-6 * max(1, abs(beta[j]))
    bp <- beta
    bm <- beta
    bp[j] <- bp[j] + h
    bm[j] <- bm[j] - h
    g[j] <- (fn(bp) - fn(bm)) / (2 * h)
  }
  g
}

delta_eval <- function(fn, beta, V, df) {
  est <- fn(beta)
  g <- num_grad(fn, beta)
  se <- sqrt(drop(t(g) %*% V %*% g))
  crit <- qt(0.975, df=df)
  c(
    estimate=est,
    se=se,
    lower=est - crit * se,
    upper=est + crit * se
  )
}

fit_config <- function(
  data,
  outcome,
  prior,
  protection,
  formula_extra="",
  case_weights=NULL,
  avg_weights=NULL,
  lower_protection=0.75,
  upper_protection=1.00,
  spec_name="spec"
) {
  d <- data
  if (any(!is.finite(d$E))) stop(paste(spec_name, "nonfinite E"))
  if (any(!is.finite(d[[protection]]))) stop(paste(spec_name, "nonfinite protection"))
  if (any(!d[[outcome]] %in% c(0,1))) stop(paste(spec_name, "invalid outcome"))
  if (any(!d[[prior]] %in% c(0,1))) stop(paste(spec_name, "invalid prior incident"))

  E_mean <- mean(d$E)
  P_mean <- mean(d[[protection]])
  d$Ec <- d$E - E_mean
  d$Pc <- d[[protection]] - P_mean

  base_terms <- c(prior, "Ec", "Pc", "Ec:Pc")
  if (nzchar(formula_extra)) base_terms <- c(base_terms, formula_extra)
  f <- as.formula(paste(outcome, "~", paste(base_terms, collapse=" + ")))

  if (is.null(case_weights)) {
    fit <- glm(f, data=d, family=binomial(link="logit"), model=TRUE, y=TRUE, na.action=na.fail)
  } else {
    if (length(case_weights) != nrow(d)) stop(paste(spec_name, "case weight length mismatch"))
    if (any(!is.finite(case_weights)) || any(case_weights <= 0)) stop(paste(spec_name, "invalid case weights"))
    fit <- suppressWarnings(glm(
      f,
      data=d,
      family=binomial(link="logit"),
      weights=case_weights,
      model=TRUE,
      y=TRUE,
      na.action=na.fail
    ))
  }

  if (!isTRUE(fit$converged)) stop(paste(spec_name, "did not converge"))
  X <- model.matrix(fit)
  if (qr(X)$rank != ncol(X)) stop(paste(spec_name, "model matrix rank failure"))

  cv <- cluster_vcov(fit, d$id)
  V <- cv$V
  beta <- coef(fit)
  if (!all(is.finite(beta))) stop(paste(spec_name, "nonfinite coefficient"))
  if (!all(is.finite(diag(V))) || any(diag(V) <= 0)) stop(paste(spec_name, "invalid variance"))

  interaction_name <- "Ec:Pc"
  if (!interaction_name %in% names(beta)) stop(paste(spec_name, "interaction term missing"))

  int_est <- unname(beta[interaction_name])
  int_se <- sqrt(V[interaction_name, interaction_name])
  int_df <- cv$G - 1
  int_crit <- qt(0.975, df=int_df)
  int_lower <- int_est - int_crit * int_se
  int_upper <- int_est + int_crit * int_se
  int_p <- 2 * pt(abs(int_est / int_se), df=int_df, lower.tail=FALSE)

  if (is.null(avg_weights)) avg_weights <- rep(1, nrow(d))
  if (length(avg_weights) != nrow(d)) stop(paste(spec_name, "average weight length mismatch"))
  avg_weights <- avg_weights / sum(avg_weights)

  make_x <- function(nd) {
    mm <- model.matrix(delete.response(terms(fit)), data=nd)
    missing_cols <- setdiff(names(beta), colnames(mm))
    if (length(missing_cols) > 0) stop(paste(spec_name, "newdata model matrix mismatch"))
    mm[, names(beta), drop=FALSE]
  }

  mean_prob <- function(b, Eval, Pval) {
    nd <- d
    nd$Ec <- Eval - E_mean
    nd$Pc <- Pval - P_mean
    mm <- make_x(nd)
    pr <- plogis(drop(mm %*% b))
    sum(avg_weights * pr)
  }

  exposure_gradient_fn <- function(Pval) {
    function(b) mean_prob(b, 1.00, Pval) - mean_prob(b, 1/3, Pval)
  }

  config_fn <- function(b) {
    exposure_gradient_fn(lower_protection)(b) - exposure_gradient_fn(upper_protection)(b)
  }

  low <- delta_eval(exposure_gradient_fn(lower_protection), beta, V, int_df)
  high <- delta_eval(exposure_gradient_fn(upper_protection), beta, V, int_df)
  conf <- delta_eval(config_fn, beta, V, int_df)

  config <- data.frame(
    spec=spec_name,
    protection_state=c("lower", "higher", "consequence"),
    protection_value=c(lower_protection, upper_protection, NA_real_),
    estimate=c(unname(low["estimate"]), unname(high["estimate"]), unname(conf["estimate"])),
    se=c(unname(low["se"]), unname(high["se"]), unname(conf["se"])),
    lower=c(unname(low["lower"]), unname(high["lower"]), unname(conf["lower"])),
    upper=c(unname(low["upper"]), unname(high["upper"]), unname(conf["upper"]))
  )

  contradiction <- isTRUE(int_lower > 0) && isTRUE(unname(conf["estimate"]) < 0)
  directional_support <- isTRUE(int_est < 0) && isTRUE(unname(conf["estimate"]) > 0)

  summary <- data.frame(
    spec=spec_name,
    estimable=TRUE,
    n=nrow(d),
    clusters=cv$G,
    events=sum(d[[outcome]]),
    interaction=int_est,
    interaction_se=int_se,
    interaction_lower=int_lower,
    interaction_upper=int_upper,
    interaction_p=int_p,
    C=unname(conf["estimate"]),
    C_se=unname(conf["se"]),
    C_lower=unname(conf["lower"]),
    C_upper=unname(conf["upper"]),
    delta_lower_protection=unname(low["estimate"]),
    delta_higher_protection=unname(high["estimate"]),
    contradiction=contradiction,
    directional_support=directional_support,
    stringsAsFactors=FALSE
  )

  list(
    fit=fit,
    V=V,
    summary=summary,
    config=config,
    data=d,
    case_weights=case_weights,
    avg_weights=avg_weights,
    E_mean=E_mean,
    P_mean=P_mean
  )
}

safe_fit <- function(...) {
  tryCatch(
    fit_config(...),
    error=function(e) list(
      fit=NULL,
      V=NULL,
      summary=data.frame(
        spec=list(...)$spec_name,
        estimable=FALSE,
        n=NA_integer_,
        clusters=NA_integer_,
        events=NA_integer_,
        interaction=NA_real_,
        interaction_se=NA_real_,
        interaction_lower=NA_real_,
        interaction_upper=NA_real_,
        interaction_p=NA_real_,
        C=NA_real_,
        C_se=NA_real_,
        C_lower=NA_real_,
        C_upper=NA_real_,
        delta_lower_protection=NA_real_,
        delta_higher_protection=NA_real_,
        contradiction=NA,
        directional_support=NA,
        error=conditionMessage(e),
        stringsAsFactors=FALSE
      ),
      config=NULL,
      error=conditionMessage(e)
    )
  )
}

primary <- fit_config(
  data=adj,
  outcome="I13_next",
  prior="I13_t",
  protection="R",
  formula_extra="factor(wave_t)",
  spec_name="primary"
)

s1 <- safe_fit(
  data=adj,
  outcome="I12_next",
  prior="I12_t",
  protection="R",
  formula_extra="factor(wave_t)",
  spec_name="S1_broader_incident"
)

s2 <- safe_fit(
  data=adj,
  outcome="I13_next",
  prior="I13_t",
  protection="R3",
  formula_extra="factor(wave_t)",
  lower_protection=2/3,
  upper_protection=1.00,
  spec_name="S2_monitoring_excluded"
)

s3 <- safe_fit(
  data=adj,
  outcome="I13_next",
  prior="I13_t",
  protection="R",
  formula_extra="factor(wave_t) + web_mode_next",
  spec_name="S3_next_wave_mode"
)

s4_data <- adj[adj$wave_t == 3, , drop=FALSE]
s4 <- safe_fit(
  data=s4_data,
  outcome="I13_next",
  prior="I13_t",
  protection="R",
  formula_extra="",
  spec_name="S4_strict_wave3_to4"
)

s5 <- safe_fit(
  data=gap2,
  outcome="I13_next",
  prior="I13_t",
  protection="R",
  formula_extra="factor(wave_t)",
  spec_name="S5_two_wave_gap"
)

s6 <- safe_fit(
  data=adj,
  outcome="I13_next",
  prior="I13_t",
  protection="R",
  formula_extra="factor(wave_t) + factor(sector)",
  spec_name="S6_sector_adjusted"
)

ret_fit <- glm(
  retain_adjacent ~ factor(wave) + I13 + E + R,
  data=elig,
  family=binomial(link="logit"),
  model=TRUE,
  y=TRUE,
  na.action=na.fail
)
if (!isTRUE(ret_fit$converged)) stop("S7 retention model did not converge")

ret_new <- data.frame(
  wave=adj$wave_t,
  I13=adj$I13_t,
  E=adj$E,
  R=adj$R
)
ret_p <- predict(ret_fit, newdata=ret_new, type="response")
if (any(!is.finite(ret_p)) || any(ret_p <= 0) || any(ret_p >= 1)) stop("S7 invalid retention probability")
ret_w <- 1 / ret_p
ret_w <- ret_w / mean(ret_w)

s7 <- safe_fit(
  data=adj,
  outcome="I13_next",
  prior="I13_t",
  protection="R",
  formula_extra="factor(wave_t)",
  case_weights=ret_w,
  avg_weights=ret_w,
  spec_name="S7_observable_retention_weighted"
)

s8 <- safe_fit(
  data=adj,
  outcome="I13_next",
  prior="I13_t",
  protection="R",
  formula_extra="factor(wave_t) + factor(region2)",
  spec_name="S8_region_adjusted"
)

sens <- list(s1, s2, s3, s4, s5, s6, s7, s8)
sens_summary <- do.call(rbind, lapply(sens, function(x) x$summary))
all_summary <- rbind(primary$summary, sens_summary)

primary_beta <- primary$summary$interaction
primary_upper <- primary$summary$interaction_upper
primary_C <- primary$summary$C
primary_low <- primary$summary$delta_lower_protection
primary_high <- primary$summary$delta_higher_protection

contradiction_flags <- sens_summary$contradiction
contradiction_flags[is.na(contradiction_flags)] <- FALSE
any_contradiction <- any(contradiction_flags)

strong_conditions <- c(
  beta_negative=primary_beta < 0,
  beta_ci_below_zero=primary_upper < 0,
  material_C=primary_C >= 0.05,
  lower_gradient_positive=primary_low > 0,
  higher_no_greater=primary_high <= primary_low,
  no_contradictory_sensitivity=!any_contradiction
)

if (all(strong_conditions)) {
  primary_decision <- "NEXT_WAVE_CONFIGURATION_SUPPORT"
} else if (
  primary_beta < 0 &&
  primary_C > 0 &&
  !any_contradiction
) {
  primary_decision <- "NEXT_WAVE_CONFIGURATION_PARTIAL"
} else {
  primary_decision <- "NEXT_WAVE_CONFIGURATION_FAIL"
}

s4_support <- isTRUE(s4$summary$estimable) && isTRUE(s4$summary$directional_support)
s5_support <- isTRUE(s5$summary$estimable) && isTRUE(s5$summary$directional_support)
s4_contra <- isTRUE(s4$summary$contradiction)
s5_contra <- isTRUE(s5$summary$contradiction)
timing_stronger_allowed <- (s4_support || s5_support) && !(s4_contra || s5_contra)

influence <- data.frame()
influence_class <- "NOT_TRIGGERED"

if (primary_decision %in% c("NEXT_WAVE_CONFIGURATION_SUPPORT", "NEXT_WAVE_CONFIGURATION_PARTIAL")) {
  b0 <- primary_beta
  ids <- unique(adj$id)
  infl_rows <- vector("list", length(ids))

  for (j in seq_along(ids)) {
    oid <- ids[j]
    dd <- adj[adj$id != oid, , drop=FALSE]
    dd$Ec <- dd$E - mean(dd$E)
    dd$Pc <- dd$R - mean(dd$R)
    f <- I13_next ~ I13_t + Ec + Pc + Ec:Pc + factor(wave_t)

    refit <- tryCatch(
      glm(f, data=dd, family=binomial(link="logit"), na.action=na.fail),
      error=function(e) NULL
    )

    if (is.null(refit) || !isTRUE(refit$converged) || !"Ec:Pc" %in% names(coef(refit))) {
      infl_rows[[j]] <- data.frame(
        id=oid,
        estimate=NA_real_,
        sign_flip=NA,
        relative_change=NA_real_,
        breach=TRUE,
        reason="nonestimable"
      )
    } else {
      bj <- unname(coef(refit)["Ec:Pc"])
      rel <- abs(bj - b0) / max(abs(b0), 1e-8)
      sign_flip <- bj > 1e-8
      breach <- isTRUE(sign_flip) || isTRUE(rel >= 0.20)
      infl_rows[[j]] <- data.frame(
        id=oid,
        estimate=bj,
        sign_flip=sign_flip,
        relative_change=rel,
        breach=breach,
        reason=ifelse(breach, "threshold", "none")
      )
    }
  }

  influence <- do.call(rbind, infl_rows)
  n_breach <- sum(influence$breach)
  if (n_breach == 0) {
    influence_class <- "PASS"
  } else if (n_breach <= 4) {
    influence_class <- "PARTIAL_TRIGGER"
  } else {
    influence_class <- "FAIL"
  }
}

retention_stats <- data.frame(
  metric=c("min","p05","median","p95","max","mean"),
  probability=c(
    min(ret_p),
    unname(quantile(ret_p, 0.05)),
    median(ret_p),
    unname(quantile(ret_p, 0.95)),
    max(ret_p),
    mean(ret_p)
  ),
  normalized_weight=c(
    min(ret_w),
    unname(quantile(ret_w, 0.05)),
    median(ret_w),
    unname(quantile(ret_w, 0.95)),
    max(ret_w),
    mean(ret_w)
  )
)

config_list <- c(list(primary$config), lapply(sens, function(x) x$config))
config_list <- config_list[!vapply(config_list, is.null, logical(1))]
all_config <- do.call(rbind, config_list)

primary_coef <- data.frame(
  term=names(coef(primary$fit)),
  estimate=coef(primary$fit),
  se=sqrt(diag(primary$V)),
  row.names=NULL
)
primary_coef$df <- primary$summary$clusters - 1
crit <- qt(0.975, df=primary_coef$df)
primary_coef$lower <- primary_coef$estimate - crit * primary_coef$se
primary_coef$upper <- primary_coef$estimate + crit * primary_coef$se
primary_coef$p_value <- 2 * pt(abs(primary_coef$estimate / primary_coef$se), df=primary_coef$df, lower.tail=FALSE)

integrity <- c(
  primary_n=nrow(adj) == 447,
  primary_orgs=length(unique(adj$id)) == 308,
  primary_events=sum(adj$I13_next) == 268,
  primary_converged=isTRUE(primary$fit$converged),
  primary_rank=qr(model.matrix(primary$fit))$rank == ncol(model.matrix(primary$fit)),
  primary_finite_coef=all(is.finite(coef(primary$fit))),
  primary_finite_vcov=all(is.finite(primary$V)),
  primary_positive_se=all(sqrt(diag(primary$V)) > 0),
  primary_config_finite=all(is.finite(primary$config$estimate)),
  s4_n=nrow(s4_data) == 89,
  s5_n=nrow(gap2) == 154,
  retention_model_converged=isTRUE(ret_fit$converged),
  retention_prob_valid=all(ret_p > 0 & ret_p < 1),
  retention_weight_valid=all(is.finite(ret_w)) && all(ret_w > 0)
)

if (!all(integrity)) {
  failed <- names(integrity)[!integrity]
  stop(paste("configuration integrity failure", paste(failed, collapse=", ")))
}

write.csv(primary_coef, file.path(out_dir, "configuration_primary_coefficients.csv"), row.names=FALSE)
write.csv(all_config, file.path(out_dir, "configuration_contrasts.csv"), row.names=FALSE)
write.csv(all_summary, file.path(out_dir, "configuration_specification_summary.csv"), row.names=FALSE)
write.csv(retention_stats, file.path(out_dir, "configuration_retention_weight_diagnostics.csv"), row.names=FALSE)
write.csv(influence, file.path(out_dir, "configuration_influence_diagnostic.csv"), row.names=FALSE)
writeLines(paste(names(integrity), integrity, sep="="), file.path(out_dir, "configuration_integrity_checks.txt"))

writeLines(c(
  paste("primary_decision", primary_decision, sep="="),
  paste("primary_beta_ER", primary_beta, sep="="),
  paste("primary_beta_lower", primary$summary$interaction_lower, sep="="),
  paste("primary_beta_upper", primary_upper, sep="="),
  paste("primary_beta_p", primary$summary$interaction_p, sep="="),
  paste("primary_C", primary_C, sep="="),
  paste("primary_C_lower", primary$summary$C_lower, sep="="),
  paste("primary_C_upper", primary$summary$C_upper, sep="="),
  paste("primary_delta_R075", primary_low, sep="="),
  paste("primary_delta_R100", primary_high, sep="="),
  paste("any_contradictory_sensitivity", any_contradiction, sep="="),
  paste("S4_directionally_supportive", s4_support, sep="="),
  paste("S5_directionally_supportive", s5_support, sep="="),
  paste("stronger_temporal_wording_allowed", timing_stronger_allowed, sep="="),
  paste("influence_class", influence_class, sep="="),
  paste("influence_breaches", ifelse(nrow(influence) > 0, sum(influence$breach), NA), sep="=")
), file.path(out_dir, "configuration_decision.txt"))

saveRDS(
  list(
    primary=primary,
    sensitivities=sens,
    retention_model=ret_fit,
    primary_decision=primary_decision,
    timing_stronger_allowed=timing_stronger_allowed,
    influence_class=influence_class
  ),
  file.path(out_dir, "configuration_model_objects.rds"),
  version=3
)

cat("CONFIGURATION_MODELS_PASS\n")
cat("primary_decision=", primary_decision, "\n", sep="")
cat("primary_beta_ER=", primary_beta, "\n", sep="")
cat("primary_beta_CI=[", primary$summary$interaction_lower, ",", primary_upper, "]\n", sep="")
cat("primary_C=", primary_C, "\n", sep="")
cat("primary_delta_R075=", primary_low, "\n", sep="")
cat("primary_delta_R100=", primary_high, "\n", sep="")
cat("timing_stronger_allowed=", timing_stronger_allowed, "\n", sep="")
cat("influence_class=", influence_class, "\n", sep="")
cat("\nSENSITIVITIES\n")
print(sens_summary)

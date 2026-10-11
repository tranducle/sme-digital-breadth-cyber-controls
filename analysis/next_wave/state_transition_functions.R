# State-dependent protection model: estimation functions.
# Pure functions only: no file I/O, no randomness.

ST_THRESHOLD_O <- 0.05
ST_THRESHOLD_S <- 0.05
ST_SIGN_TOL <- 1e-8
ST_INFLUENCE_REL <- 0.20
ST_INFLUENCE_FAIL_COUNT <- 5L
ST_CI_LEVEL <- 0.95
ST_VERDICT_SUPPORT <- "STATE_DEPENDENT_PROTECTION_SUPPORT"
ST_VERDICT_PARTIAL <- "STATE_DEPENDENT_PROTECTION_PARTIAL"
ST_VERDICT_FAIL <- "STATE_DEPENDENT_PROTECTION_FAIL"
ST_FOCAL_TERMS <- c("It", "Rc", "Ec", "It:Rc", "It:Ec")
ST_SEPARATION_PATTERN <- "fitted probabilities numerically 0 or 1"

# Build one analysis sample. Centering is within the analysis sample; beta_R, beta_IR and
# every standardized contrast are invariant to the centering constant because the model
# carries both It:Rc and It:Ec (the shift is absorbed by alpha and beta_I).
st_build_sample <- function(df, incident, protection, timing, extras = character(0)) {
  stopifnot(incident %in% c("I13", "I12"), protection %in% c("R", "R3"),
            timing %in% c("adjacent", "none", "gap2"))
  d <- data.frame(
    id = df$id,
    y = as.numeric(df[[paste0(incident, "_next")]]),
    It = as.numeric(df[[paste0(incident, "_t")]]),
    Rraw = as.numeric(df[[protection]]),
    E = as.numeric(df$E),
    wave_t = as.integer(df$wave_t)
  )
  r_center <- mean(d$Rraw)
  e_center <- mean(d$E)
  d$Rc <- d$Rraw - r_center
  d$Ec <- d$E - e_center
  terms <- character(0)
  if (timing == "adjacent") {
    d$T23 <- as.numeric(d$wave_t == 2)
    d$T34 <- as.numeric(d$wave_t == 3)
    terms <- c(terms, "T23", "T34")
  }
  if (timing == "gap2") {
    d$T24 <- as.numeric(d$wave_t == 2)
    terms <- c(terms, "T24")
  }
  if ("web_mode_next" %in% extras) {
    d$web_mode_next <- as.numeric(df$web_mode_next)
    terms <- c(terms, "web_mode_next")
  }
  if ("sector" %in% extras) {
    d$sector_f <- factor(df$sector, levels = sort(unique(df$sector)))
    terms <- c(terms, "sector_f")
  }
  if ("region2" %in% extras) {
    d$region_f <- factor(df$region2, levels = sort(unique(df$region2)))
    terms <- c(terms, "region_f")
  }
  rhs <- paste(c("It", "Rc", "Ec", "It:Rc", "It:Ec", terms), collapse = " + ")
  list(data = d, formula = as.formula(paste("y ~", rhs)),
       r_center = r_center, e_center = e_center)
}

# Fit the binomial logit glm, capturing (not suppressing) every warning.
st_fit_glm <- function(formula, data, weights = NULL) {
  seen <- character(0)
  data$.c5w <- if (is.null(weights)) rep(1, nrow(data)) else as.numeric(weights)
  fit <- tryCatch(
    withCallingHandlers(
      glm(formula, data = data, family = binomial(link = "logit"), weights = .c5w,
          model = TRUE, y = TRUE, na.action = na.fail),
      warning = function(w) {
        seen <<- c(seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) structure(list(message = conditionMessage(e)), class = "st_fit_error")
  )
  list(fit = fit, warnings = unique(seen),
       separation = any(grepl(ST_SEPARATION_PATTERN, seen, fixed = TRUE)))
}

# Organization-clustered score sandwich with the CR1 multiplier G/(G-1) * (N-1)/(N-K).
# bread %*% crossprod(S) %*% bread is computed as crossprod(S %*% bread): the same algebra, but
# every diagonal entry is a sum of squares. In the first order the diagonal of an ill-conditioned
# fit can come out negative by round-off alone (synthetic S6 under quasi-separated sector levels,
# information condition 1.0e10: intercept variance -1.89 in numpy, +1.03 in R), which would make
# estimability depend on the arithmetic of one implementation. The information matrix is inverted
# by solve() with its default rule (LAPACK dgesv; an exactly singular system or a reciprocal
# condition number below machine epsilon is an error); no stricter rank rule is added (v2 lines
# 161 and 177 define estimability by model-matrix rank and finite covariance values).
st_cluster_vcov <- function(model, cluster) {
  X <- model.matrix(model)
  y <- model.response(model.frame(model))
  mu <- fitted(model)
  pw <- model$prior.weights
  if (is.null(pw)) pw <- rep(1, length(y))
  if (length(cluster) != length(y)) stop("cluster length mismatch")
  score_cluster <- rowsum(X * as.numeric(pw * (y - mu)), group = cluster, reorder = FALSE)
  info <- crossprod(X, X * as.numeric(pw * mu * (1 - mu)))
  bread <- solve(info)
  N <- nrow(X)
  K <- ncol(X)
  G <- length(unique(cluster))
  if (G < 2) stop("fewer than two clusters")
  cr1 <- (G / (G - 1)) * ((N - 1) / (N - K))
  V <- cr1 * crossprod(score_cluster %*% bread)
  dimnames(V) <- list(colnames(X), colnames(X))
  if (!all(is.finite(V))) stop("nonfinite cluster covariance")
  list(V = V, G = G, N = N, K = K, cr1 = cr1)
}

# Estimability checks used for every specification and every leave-one-out refit.
st_check_fit <- function(fitres, cluster, need_vcov = TRUE) {
  fit <- fitres$fit
  if (inherits(fit, "st_fit_error")) {
    return(list(ok = FALSE, reasons = paste("glm error:", fit$message), vc = NULL))
  }
  reasons <- character(0)
  if (!isTRUE(fit$converged)) reasons <- c(reasons, "glm did not converge")
  X <- model.matrix(fit)
  if (qr(X)$rank != ncol(X)) reasons <- c(reasons, "model matrix rank deficient")
  b <- coef(fit)
  if (any(!is.finite(b))) reasons <- c(reasons, "non-finite or aliased coefficient")
  absent <- setdiff(ST_FOCAL_TERMS, names(b))
  if (length(absent)) reasons <- c(reasons, paste("missing focal term:", paste(absent, collapse = ",")))
  if (!length(reasons)) {
    eta_gap <- max(abs(drop(X %*% b) - fit$linear.predictors))
    if (!is.finite(eta_gap) || eta_gap > 1e-8) reasons <- c(reasons, "linear predictor not reproduced")
  }
  vc <- NULL
  if (!length(reasons) && need_vcov) {
    vc <- tryCatch(st_cluster_vcov(fit, cluster), error = function(e) conditionMessage(e))
    if (is.character(vc)) {
      reasons <- c(reasons, paste("covariance failure:", vc))
      vc <- NULL
    } else if (any(!is.finite(sqrt(diag(vc$V)))) || any(diag(vc$V) <= 0)) {
      reasons <- c(reasons, "non-positive or non-finite standard error")
      vc <- NULL
    }
  }
  list(ok = !length(reasons), reasons = reasons, vc = vc)
}

st_t_crit <- function(G) qt(1 - (1 - ST_CI_LEVEL) / 2, df = G - 1)

st_coef_table <- function(fit, vc) {
  b <- coef(fit)
  se <- sqrt(diag(vc$V))[names(b)]
  tc <- st_t_crit(vc$G)
  stat <- b / se
  data.frame(term = names(b), estimate = unname(b), se = unname(se),
             lcl = unname(b - tc * se), ucl = unname(b + tc * se),
             t_stat = unname(stat), p_value = unname(2 * pt(-abs(stat), df = vc$G - 1)),
             df = vc$G - 1, stringsAsFactors = FALSE)
}

# Counterfactual design rows: set It and raw protection, keep E, timing and factor columns.
st_counterfactual <- function(X, set_it, r_value, r_center) {
  Xc <- X
  rc <- r_value - r_center
  Xc[, "It"] <- set_it
  Xc[, "Rc"] <- rc
  Xc[, "It:Rc"] <- set_it * rc
  Xc[, "It:Ec"] <- set_it * Xc[, "Ec"]
  Xc
}

st_std_prob <- function(Xc, b, w) {
  p <- plogis(drop(Xc %*% b))
  wn <- w / sum(w)
  list(value = sum(wn * p), grad = colSums(Xc * as.numeric(wn * p * (1 - p))))
}

# Standardized contrast Pr(R = r_low) - Pr(R = r_high) over `rows`, with delta-method gradient.
st_contrast <- function(fit, rows, set_it, r_low, r_high, r_center, avg_w) {
  X <- model.matrix(fit)[rows, , drop = FALSE]
  b <- coef(fit)
  w <- avg_w[rows]
  lo <- st_std_prob(st_counterfactual(X, set_it, r_low, r_center), b, w)
  hi <- st_std_prob(st_counterfactual(X, set_it, r_high, r_center), b, w)
  list(estimate = lo$value - hi$value, grad = lo$grad - hi$grad,
       p_low = lo$value, p_high = hi$value, n_avg = length(rows))
}

st_contrast_row <- function(name, ct, V, G) {
  se <- sqrt(drop(t(ct$grad) %*% V %*% ct$grad))
  if (!is.finite(se) || se <= 0) stop(paste("non-finite contrast SE for", name))
  tc <- st_t_crit(G)
  data.frame(contrast = name, estimate = ct$estimate, se = se,
             lcl = ct$estimate - tc * se, ucl = ct$estimate + tc * se,
             prob_low_state = ct$p_low, prob_high_state = ct$p_high,
             n_averaged = ct$n_avg, stringsAsFactors = FALSE)
}

# All five contrasts (O, P, O_common, P_common, S) for one fitted specification.
st_all_contrasts <- function(fit, vc, it_vec, r_low, r_high, r_center, avg_w) {
  rows0 <- which(it_vec == 0)
  rows1 <- which(it_vec == 1)
  rows_all <- seq_along(it_vec)
  if (!length(rows0) || !length(rows1)) stop("an incident state has no rows to average over")
  O <- st_contrast(fit, rows0, 0, r_low, r_high, r_center, avg_w)
  P <- st_contrast(fit, rows1, 1, r_low, r_high, r_center, avg_w)
  Oc <- st_contrast(fit, rows_all, 0, r_low, r_high, r_center, avg_w)
  Pc <- st_contrast(fit, rows_all, 1, r_low, r_high, r_center, avg_w)
  S <- list(estimate = Oc$estimate - Pc$estimate, grad = Oc$grad - Pc$grad,
            p_low = NA_real_, p_high = NA_real_, n_avg = length(rows_all))
  rbind(st_contrast_row("O", O, vc$V, vc$G), st_contrast_row("P", P, vc$V, vc$G),
        st_contrast_row("O_common", Oc, vc$V, vc$G), st_contrast_row("P_common", Pc, vc$V, vc$G),
        st_contrast_row("S", S, vc$V, vc$G))
}

# Contradiction rule (applied uniformly over S1 to S8):
# incident-free protection slope statistically positive AND onset contrast negative.
st_is_contradictory <- function(betaR_lcl, O) {
  isTRUE(betaR_lcl > 0) && isTRUE(O < 0)
}

# The nine SUPPORT conditions of the decision rule, as a named logical vector.
st_support_conditions <- function(bR, bR_ucl, O, O_lcl, bIR, bIR_lcl, bIR_ucl,
                                  S, S_lcl, S_ucl, any_contradiction, tol = ST_SIGN_TOL) {
  c(C1_betaR_negative = bR < 0,
    C2_betaR_upper_below_0 = bR_ucl < 0,
    C3_O_at_least_threshold = O >= ST_THRESHOLD_O,
    C4_O_lower_above_0 = O_lcl > 0,
    ST_betaIR_ci_excludes_0 = (bIR_lcl > 0) || (bIR_ucl < 0),
    C6_absS_at_least_threshold = abs(S) >= ST_THRESHOLD_S,
    C7_S_ci_excludes_0 = (S_lcl > 0) || (S_ucl < 0),
    C8_betaIR_S_same_nonzero_sign = abs(bIR) > tol && abs(S) > tol && sign(bIR) == sign(S),
    C9_no_contradictory_sensitivity = !any_contradiction)
}

st_classify <- function(conditions, bR, O, integrity_ok) {
  if (!isTRUE(integrity_ok) || anyNA(conditions)) return(ST_VERDICT_FAIL)
  if (all(conditions)) return(ST_VERDICT_SUPPORT)
  if (bR < 0 && O > 0 && isTRUE(conditions[["C9_no_contradictory_sensitivity"]])) {
    return(ST_VERDICT_PARTIAL)
  }
  ST_VERDICT_FAIL
}

st_rel_change <- function(full, loo, tol) {
  if (abs(full) <= tol) return(if (abs(loo - full) > tol) Inf else 0)
  abs(loo - full) / abs(full)
}

# Influence breach rule for one leave-one-organization-out refit.
st_influence_breach <- function(full_bR, full_bIR, loo_bR, loo_bIR, estimable,
                                tol = ST_SIGN_TOL, rel = ST_INFLUENCE_REL) {
  if (!isTRUE(estimable)) {
    return(c(betaR_neg_to_pos = NA, betaIR_sign_change = NA, betaR_rel_ge = NA,
             betaIR_rel_ge = NA, non_estimable = TRUE, breach = TRUE))
  }
  f <- c(betaR_neg_to_pos = full_bR < 0 && loo_bR > tol,
         betaIR_sign_change = (full_bIR > tol && loo_bIR < -tol) || (full_bIR < -tol && loo_bIR > tol),
         betaR_rel_ge = st_rel_change(full_bR, loo_bR, tol) >= rel,
         betaIR_rel_ge = st_rel_change(full_bIR, loo_bIR, tol) >= rel,
         non_estimable = FALSE)
  c(f, breach = any(f))
}

st_influence_class <- function(n_breach) {
  if (n_breach == 0) return("PASS")
  if (n_breach < ST_INFLUENCE_FAIL_COUNT) return("PARTIAL_TRIGGER")
  "FAIL"
}

# Leave-one-organization-out refits of the primary specification.
st_influence <- function(sample, full_bR, full_bIR) {
  ids <- unique(sample$data$id)
  rows <- lapply(ids, function(g) {
    d <- sample$data[sample$data$id != g, , drop = FALSE]
    fr <- st_fit_glm(sample$formula, d)
    chk <- st_check_fit(fr, d$id, need_vcov = FALSE)
    loo_bR <- if (chk$ok) unname(coef(fr$fit)["Rc"]) else NA_real_
    loo_bIR <- if (chk$ok) unname(coef(fr$fit)["It:Rc"]) else NA_real_
    br <- st_influence_breach(full_bR, full_bIR, loo_bR, loo_bIR, chk$ok)
    data.frame(dropped_id = g, n_rows_dropped = sum(sample$data$id == g), estimable = chk$ok,
               betaR = loo_bR, betaIR = loo_bIR,
               betaR_rel_change = if (chk$ok) st_rel_change(full_bR, loo_bR, ST_SIGN_TOL) else NA_real_,
               betaIR_rel_change = if (chk$ok) st_rel_change(full_bIR, loo_bIR, ST_SIGN_TOL) else NA_real_,
               separation_warning = fr$separation, t(br), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

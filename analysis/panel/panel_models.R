epsilon <- 1e-10
cov_epsilon <- 1e-12
machine_epsilon <- 2.220446049250313e-16

if (as.character(getRversion()) != "4.6.1") {
  warning("Results were produced with R 4.6.1; running ", as.character(getRversion()))
}

root <- "."
out_dir <- "results/panel"
raw_csv <- file.path(out_dir, "csls_raw_selected_columns.csv")
source_path <- file.path(root, "data/raw/UKDA-9356-stata/stata/stata13/csls_waves_1-4.dta")

expected_source_sha256 <- "ef4632c26fb1a02c3278bafb707ee11c1a8acb26c7a33d7c2dedcc272c71a306"
expected_raw_csv_sha256 <- "d236650ae3576470b82db13586975f8f28d1d67175864d660013de7b294faa87"

sha256 <- function(path) {
  out <- system2("shasum", c("-a", "256", path), stdout=TRUE, stderr=TRUE)
  if (length(out) < 1) stop("sha256 command failed")
  strsplit(out[1], "[[:space:]]+")[[1]][1]
}

if (sha256(source_path) != expected_source_sha256) stop("source dataset hash drift")
if (sha256(raw_csv) != expected_raw_csv_sha256) stop("raw bridge csv hash drift")

col_classes <- c(
  source_row="integer", wave="numeric", qtypdum="numeric", LongId="character",
  onlinea="numeric", onlineb="numeric", devices="numeric",
  rulesb="numeric", rulesf="numeric", rulesh="numeric", rulesi="numeric",
  qsize="numeric", modetype="numeric"
)

raw <- read.csv(
  raw_csv,
  stringsAsFactors=FALSE,
  check.names=FALSE,
  colClasses=col_classes,
  na.strings="NA"
)

required_cols <- names(col_classes)
if (!identical(names(raw), required_cols)) stop("raw csv columns drift")
if (nrow(raw) != 4876) stop("raw row count drift")
if (anyDuplicated(raw$source_row)) stop("source row identifier duplicated")

raw$id <- trimws(raw$LongId)
binary_items <- c("onlinea","onlineb","devices","rulesb","rulesf","rulesh","rulesi")
valid_binary <- function(x) !is.na(x) & x %in% c(1,2)

eligible <- raw$qtypdum == 1 &
  raw$qsize == 2 &
  !is.na(raw$id) &
  nzchar(raw$id)

for (v in binary_items) eligible <- eligible & valid_binary(raw[[v]])

d <- raw[eligible, c("source_row","wave","id","modetype",binary_items), drop=FALSE]

for (v in binary_items) d[[paste0(v,"_bin")]] <- as.numeric(d[[v]] == 1)

d$E_score <- d$onlinea_bin + d$onlineb_bin + d$devices_bin
d$E <- d$E_score / 3

d$R_score <- d$rulesb_bin + d$rulesf_bin + d$rulesh_bin + d$rulesi_bin
d$R <- d$R_score / 4

d$R3_score <- d$rulesf_bin + d$rulesh_bin + d$rulesi_bin
d$R3 <- d$R3_score / 3

d$web_mode <- as.numeric(d$modetype == 2)

org_n <- table(d$id)
keep_ids <- names(org_n)[org_n >= 2]
panel <- d[d$id %in% keep_ids, , drop=FALSE]
panel <- panel[order(panel$id, panel$wave, panel$source_row), , drop=FALSE]
row.names(panel) <- NULL

if (anyDuplicated(paste(panel$id, panel$wave, sep="|"))) stop("duplicate organization wave in primary panel")
if (nrow(panel) != 783) stop("panel row count drift")
if (length(unique(panel$id)) != 318) stop("organization count drift")

org_var_count <- function(x, id) {
  z <- tapply(x, id, function(v) max(v) - min(v) > epsilon)
  sum(z)
}

n_e_var <- org_var_count(panel$E, panel$id)
n_r_var <- org_var_count(panel$R, panel$id)
both_var_by_org <- tapply(seq_len(nrow(panel)), panel$id, function(ii) {
  (max(panel$E[ii]) - min(panel$E[ii]) > epsilon) &&
    (max(panel$R[ii]) - min(panel$R[ii]) > epsilon)
})
n_both_var <- sum(both_var_by_org)

if (n_e_var != 188) stop("exposure variation count drift")
if (n_r_var != 128) stop("response variation count drift")
if (n_both_var != 89) stop("joint variation count drift")

write.csv(panel, file.path(out_dir, "analysis_ready_panel.csv"), row.names=FALSE)

rank_info <- function(X) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (ncol(X) == 0 || nrow(X) == 0) {
    return(list(ok=FALSE, rank=0, tol=NA_real_, smax=NA_real_, reason="empty_matrix"))
  }
  if (any(!is.finite(X))) {
    return(list(ok=FALSE, rank=NA_integer_, tol=NA_real_, smax=NA_real_, reason="nonfinite_matrix"))
  }
  sv <- tryCatch(svd(X, nu=0, nv=0)$d, error=function(e) NULL)
  if (is.null(sv) || length(sv) == 0 || any(!is.finite(sv))) {
    return(list(ok=FALSE, rank=NA_integer_, tol=NA_real_, smax=NA_real_, reason="svd_failure"))
  }
  smax <- max(sv)
  tol <- max(nrow(X), ncol(X)) * machine_epsilon * smax
  rk <- sum(sv > tol)
  list(ok=TRUE, rank=rk, tol=tol, smax=smax, reason="ok")
}

demean_vec <- function(x, id) {
  x - ave(x, id, FUN=mean)
}

cluster_fit <- function(X, y, cluster, coef_names, model_name) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  cluster <- as.character(cluster)
  N <- nrow(X)
  K <- ncol(X)
  G <- length(unique(cluster))
  ri <- rank_info(X)

  base_diag <- list(
    model=model_name, N=N, K=K, G=G,
    rank=ri$rank, rank_tol=ri$tol, rank_reason=ri$reason,
    cr1_multiplier=NA_real_, numerical_status="not_estimated"
  )

  fail <- function(reason) {
    base_diag$numerical_status <<- reason
    list(ok=FALSE, reason=reason, diagnostics=base_diag)
  }

  if (!ri$ok) return(fail(paste0("rank_check_",ri$reason)))
  if (ri$rank != K) return(fail("rank_deficient"))
  if (G <= 1) return(fail("invalid_cluster_count"))
  if (N <= K) return(fail("invalid_residual_dimension"))
  if (any(!is.finite(y))) return(fail("nonfinite_outcome"))

  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  beta <- tryCatch(as.numeric(solve(XtX, Xty)), error=function(e) NULL)
  if (is.null(beta) || length(beta) != K || any(!is.finite(beta))) return(fail("coefficient_failure"))

  u <- as.numeric(y - X %*% beta)
  if (any(!is.finite(u))) return(fail("residual_failure"))

  M <- matrix(0, K, K)
  for (g in unique(cluster)) {
    ii <- which(cluster == g)
    sg <- tryCatch(crossprod(X[ii,,drop=FALSE], u[ii]), error=function(e) NULL)
    if (is.null(sg) || any(!is.finite(sg))) return(fail("cluster_score_failure"))
    M <- M + sg %*% t(sg)
  }

  B <- tryCatch(solve(XtX), error=function(e) NULL)
  if (is.null(B) || any(!is.finite(B))) return(fail("bread_failure"))

  denom <- N - K
  if (denom <= 0 || G - 1 <= 0) return(fail("invalid_cr1_denominator"))
  cr1 <- (G/(G-1)) * ((N-1)/denom)
  if (!is.finite(cr1) || cr1 <= 0) return(fail("invalid_cr1_multiplier"))

  V <- cr1 * (B %*% M %*% B)
  if (any(!is.finite(V))) return(fail("nonfinite_covariance"))

  vd <- diag(V)
  if (any(vd < -cov_epsilon)) return(fail("negative_covariance_diagonal"))
  tiny_neg <- vd < 0 & vd >= -cov_epsilon
  if (any(tiny_neg)) {
    vd[tiny_neg] <- 0
    diag(V) <- vd
  }

  se <- sqrt(diag(V))
  if (any(!is.finite(se))) return(fail("standard_error_failure"))

  crit <- qt(0.975, df=G-1)
  if (!is.finite(crit)) return(fail("critical_value_failure"))

  lower <- beta - crit * se
  upper <- beta + crit * se
  if (any(!is.finite(lower)) || any(!is.finite(upper))) return(fail("confidence_interval_failure"))

  names(beta) <- coef_names
  names(se) <- coef_names
  names(lower) <- coef_names
  names(upper) <- coef_names
  dimnames(V) <- list(coef_names, coef_names)

  base_diag$cr1_multiplier <- cr1
  base_diag$numerical_status <- "PASS"

  list(
    ok=TRUE, reason="ok", beta=beta, V=V, se=se,
    lower=lower, upper=upper, diagnostics=base_diag
  )
}

primary_matrix <- function(dat) {
  id <- dat$id
  w2 <- as.numeric(dat$wave == 2)
  w3 <- as.numeric(dat$wave == 3)
  w4 <- as.numeric(dat$wave == 4)
  X <- cbind(
    E=demean_vec(dat$E, id),
    Wave2=demean_vec(w2, id),
    Wave3=demean_vec(w3, id),
    Wave4=demean_vec(w4, id)
  )
  y <- demean_vec(dat$R, id)
  list(X=X, y=y)
}

diagnostic_rows <- list()
add_diag <- function(fit) {
  diagnostic_rows[[length(diagnostic_rows)+1]] <<- as.data.frame(fit$diagnostics, stringsAsFactors=FALSE)
}

sample_counts <- data.frame(
  metric=c(
    "rows","organizations","organizations_exposure_variation",
    "organizations_response_variation","organizations_both_variation"
  ),
  value=c(nrow(panel),length(unique(panel$id)),n_e_var,n_r_var,n_both_var)
)
write.csv(sample_counts, file.path(out_dir, "sample_and_variation_counts.csv"), row.names=FALSE)

primary_fail_floor <- n_e_var < 100 || n_both_var < 50
pm <- primary_matrix(panel)
primary_fit <- cluster_fit(
  pm$X, pm$y, panel$id,
  c("E","Wave2","Wave3","Wave4"),
  "primary"
)
add_diag(primary_fit)

if (primary_fail_floor || !primary_fit$ok) {
  primary_status <- "FAIL"
} else {
  b <- primary_fit$beta["E"]
  lo <- primary_fit$lower["E"]
  if (b <= epsilon) {
    primary_status <- "FAIL"
  } else if (lo <= epsilon) {
    primary_status <- "PARTIAL"
  } else {
    primary_status <- "SUPPORT"
  }
}

primary_table <- if (primary_fit$ok) {
  data.frame(
    term=names(primary_fit$beta),
    estimate=as.numeric(primary_fit$beta),
    se=as.numeric(primary_fit$se),
    lower=as.numeric(primary_fit$lower),
    upper=as.numeric(primary_fit$upper)
  )
} else {
  data.frame(term="E",estimate=NA_real_,se=NA_real_,lower=NA_real_,upper=NA_real_)
}
write.csv(primary_table, file.path(out_dir, "primary_coefficients.csv"), row.names=FALSE)
if (primary_fit$ok) write.csv(primary_fit$V, file.path(out_dir, "primary_CR1_covariance.csv"), row.names=TRUE)

make_adjacent_pairs <- function(dat) {
  pieces <- list()
  for (ident in unique(dat$id)) {
    g <- dat[dat$id == ident, , drop=FALSE]
    g <- g[order(g$wave), , drop=FALSE]
    if (nrow(g) < 2) next
    for (j in 2:nrow(g)) {
      if (g$wave[j] - g$wave[j-1] != 1) next
      pieces[[length(pieces)+1]] <- data.frame(
        id=ident,
        wave0=g$wave[j-1],
        wave1=g$wave[j],
        dE=g$E[j]-g$E[j-1],
        dR=g$R[j]-g$R[j-1],
        R0_score=g$R_score[j-1],
        R1_score=g$R_score[j]
      )
    }
  }
  if (length(pieces) == 0) return(data.frame())
  do.call(rbind, pieces)
}

pairs <- make_adjacent_pairs(panel)
write.csv(pairs, file.path(out_dir, "adjacent_pairs.csv"), row.names=FALSE)

org_nonzero_de <- length(unique(pairs$id[abs(pairs$dE) > epsilon]))
org_both_delta <- length(unique(pairs$id[abs(pairs$dE) > epsilon & abs(pairs$dR) > epsilon]))
s1_floor_fail <- org_nonzero_de < 100 || org_both_delta < 50

if (nrow(pairs) > 0) {
  X1 <- cbind(
    Intercept=1,
    dE=pairs$dE,
    Transition23=as.numeric(pairs$wave0 == 2 & pairs$wave1 == 3),
    Transition34=as.numeric(pairs$wave0 == 3 & pairs$wave1 == 4)
  )
  s1_fit <- cluster_fit(
    X1, pairs$dR, pairs$id,
    c("Intercept","dE","Transition23","Transition34"),
    "S1_first_difference"
  )
} else {
  s1_fit <- list(
    ok=FALSE,
    reason="no_adjacent_pairs",
    diagnostics=list(
      model="S1_first_difference",N=0,K=4,G=0,rank=NA_integer_,
      rank_tol=NA_real_,rank_reason="empty_matrix",cr1_multiplier=NA_real_,
      numerical_status="no_adjacent_pairs"
    )
  )
}
add_diag(s1_fit)

if (s1_floor_fail || !s1_fit$ok) {
  s1_status <- "FAIL"
} else {
  b <- s1_fit$beta["dE"]
  up <- s1_fit$upper["dE"]
  if (up < -epsilon) {
    s1_status <- "FAIL"
  } else if (b < -epsilon) {
    s1_status <- "PARTIAL_TRIGGER"
  } else {
    s1_status <- "PASS"
  }
}

s1_table <- if (s1_fit$ok) {
  data.frame(
    term=names(s1_fit$beta),
    estimate=as.numeric(s1_fit$beta),
    se=as.numeric(s1_fit$se),
    lower=as.numeric(s1_fit$lower),
    upper=as.numeric(s1_fit$upper)
  )
} else {
  data.frame(term="dE",estimate=NA_real_,se=NA_real_,lower=NA_real_,upper=NA_real_)
}
write.csv(s1_table, file.path(out_dir, "S1_first_difference.csv"), row.names=FALSE)

r3_var_count <- org_var_count(panel$R3, panel$id)
r3_var_prop <- r3_var_count / length(unique(panel$id))

if (r3_var_prop < 0.10) {
  s2_status <- "NON_IDENTIFYING"
  s2_fit <- list(
    ok=FALSE,
    reason="below_variation_floor",
    diagnostics=list(
      model="S2_monitoring_excluded",N=nrow(panel),K=4,G=length(unique(panel$id)),
      rank=rank_info(pm$X)$rank,rank_tol=rank_info(pm$X)$tol,
      rank_reason="not_fit_below_variation_floor",cr1_multiplier=NA_real_,
      numerical_status="NON_IDENTIFYING"
    )
  )
} else {
  y2 <- demean_vec(panel$R3, panel$id)
  s2_fit <- cluster_fit(
    pm$X, y2, panel$id,
    c("E","Wave2","Wave3","Wave4"),
    "S2_monitoring_excluded"
  )
  if (!s2_fit$ok) {
    s2_status <- "FAIL"
  } else {
    b <- s2_fit$beta["E"]
    up <- s2_fit$upper["E"]
    if (up < -epsilon) {
      s2_status <- "FAIL"
    } else if (b < -epsilon) {
      s2_status <- "PARTIAL_TRIGGER"
    } else {
      s2_status <- "PASS"
    }
  }
}
add_diag(s2_fit)

s2_table <- if (s2_fit$ok) {
  data.frame(
    term=names(s2_fit$beta),
    estimate=as.numeric(s2_fit$beta),
    se=as.numeric(s2_fit$se),
    lower=as.numeric(s2_fit$lower),
    upper=as.numeric(s2_fit$upper)
  )
} else {
  data.frame(term="E",estimate=NA_real_,se=NA_real_,lower=NA_real_,upper=NA_real_)
}
write.csv(s2_table, file.path(out_dir, "S2_monitoring_excluded.csv"), row.names=FALSE)

id <- panel$id
w2 <- as.numeric(panel$wave == 2)
w3 <- as.numeric(panel$wave == 3)
w4 <- as.numeric(panel$wave == 4)
X3 <- cbind(
  CloudServer=demean_vec(panel$onlinea_bin,id),
  PhysicalServer=demean_vec(panel$onlineb_bin,id),
  PersonalDeviceAccess=demean_vec(panel$devices_bin,id),
  Wave2=demean_vec(w2,id),
  Wave3=demean_vec(w3,id),
  Wave4=demean_vec(w4,id)
)
y3 <- demean_vec(panel$R,id)
s3_fit <- cluster_fit(
  X3,y3,id,
  c("CloudServer","PhysicalServer","PersonalDeviceAccess","Wave2","Wave3","Wave4"),
  "S3_domain_decomposition"
)
add_diag(s3_fit)

if (!s3_fit$ok) {
  s3_status <- "FAIL"
} else {
  neg_ci <- s3_fit$upper[c("CloudServer","PhysicalServer","PersonalDeviceAccess")] < -epsilon
  n_neg_ci <- sum(neg_ci)
  if (n_neg_ci >= 2) {
    s3_status <- "FAIL"
  } else if (n_neg_ci == 1) {
    s3_status <- "PARTIAL_TRIGGER"
  } else {
    s3_status <- "PASS"
  }
}

s3_table <- if (s3_fit$ok) {
  data.frame(
    term=names(s3_fit$beta),
    estimate=as.numeric(s3_fit$beta),
    se=as.numeric(s3_fit$se),
    lower=as.numeric(s3_fit$lower),
    upper=as.numeric(s3_fit$upper)
  )
} else {
  data.frame(term="NON_ESTIMABLE",estimate=NA_real_,se=NA_real_,lower=NA_real_,upper=NA_real_)
}
write.csv(s3_table, file.path(out_dir, "S3_domain_decomposition.csv"), row.names=FALSE)

m <- demean_vec(panel$web_mode, panel$id)
primary_rank <- rank_info(pm$X)$rank
if (max(abs(m)) <= epsilon) {
  s4_status <- "NON_IDENTIFYING"
  s4_fit <- list(
    ok=FALSE,
    reason="mode_no_within_variation",
    diagnostics=list(
      model="S4_mode_adjustment",N=nrow(panel),K=5,G=length(unique(panel$id)),
      rank=primary_rank,rank_tol=rank_info(pm$X)$tol,
      rank_reason="mode_no_within_variation",cr1_multiplier=NA_real_,
      numerical_status="NON_IDENTIFYING"
    )
  )
} else {
  X4 <- cbind(pm$X, WebMode=m)
  r4 <- rank_info(X4)
  if (!r4$ok) {
    s4_status <- "FAIL"
    s4_fit <- list(
      ok=FALSE,reason="rank_check_failure",
      diagnostics=list(
        model="S4_mode_adjustment",N=nrow(panel),K=5,G=length(unique(panel$id)),
        rank=r4$rank,rank_tol=r4$tol,rank_reason=r4$reason,
        cr1_multiplier=NA_real_,numerical_status="rank_check_failure"
      )
    )
  } else if (r4$rank == primary_rank) {
    s4_status <- "NON_IDENTIFYING"
    s4_fit <- list(
      ok=FALSE,reason="mode_redundant",
      diagnostics=list(
        model="S4_mode_adjustment",N=nrow(panel),K=5,G=length(unique(panel$id)),
        rank=r4$rank,rank_tol=r4$tol,rank_reason="mode_redundant",
        cr1_multiplier=NA_real_,numerical_status="NON_IDENTIFYING"
      )
    )
  } else if (r4$rank == 5) {
    s4_fit <- cluster_fit(
      X4,pm$y,panel$id,
      c("E","Wave2","Wave3","Wave4","WebMode"),
      "S4_mode_adjustment"
    )
    add_diag(s4_fit)
    if (!s4_fit$ok) {
      s4_status <- "FAIL"
    } else {
      b <- s4_fit$beta["E"]
      up <- s4_fit$upper["E"]
      if (up < -epsilon) {
        s4_status <- "FAIL"
      } else if (b < -epsilon) {
        s4_status <- "PARTIAL_TRIGGER"
      } else {
        s4_status <- "PASS"
      }
    }
  } else {
    s4_status <- "FAIL"
    s4_fit <- list(
      ok=FALSE,reason="unexpected_augmented_rank",
      diagnostics=list(
        model="S4_mode_adjustment",N=nrow(panel),K=5,G=length(unique(panel$id)),
        rank=r4$rank,rank_tol=r4$tol,rank_reason="unexpected_augmented_rank",
        cr1_multiplier=NA_real_,numerical_status="FAIL"
      )
    )
  }
}
if (!any(vapply(diagnostic_rows, function(z) identical(as.character(z$model[1]),"S4_mode_adjustment"), logical(1)))) {
  add_diag(s4_fit)
}

s4_table <- if (isTRUE(s4_fit$ok)) {
  data.frame(
    term=names(s4_fit$beta),
    estimate=as.numeric(s4_fit$beta),
    se=as.numeric(s4_fit$se),
    lower=as.numeric(s4_fit$lower),
    upper=as.numeric(s4_fit$upper)
  )
} else {
  data.frame(term=s4_status,estimate=NA_real_,se=NA_real_,lower=NA_real_,upper=NA_real_)
}
write.csv(s4_table, file.path(out_dir, "S4_mode_adjustment.csv"), row.names=FALSE)

coef_primary_only <- function(dat) {
  mm <- primary_matrix(dat)
  ri <- rank_info(mm$X)
  if (!ri$ok || ri$rank != ncol(mm$X)) {
    return(list(ok=FALSE,beta=NA_real_,rank=ri$rank,rank_tol=ri$tol,reason="rank_failure"))
  }
  b <- tryCatch(as.numeric(solve(crossprod(mm$X), crossprod(mm$X,mm$y))), error=function(e) NULL)
  if (is.null(b) || length(b) < 1 || any(!is.finite(b))) {
    return(list(ok=FALSE,beta=NA_real_,rank=ri$rank,rank_tol=ri$tol,reason="coefficient_failure"))
  }
  list(ok=TRUE,beta=b[1],rank=ri$rank,rank_tol=ri$tol,reason="ok")
}

influence_rows <- list()
if (primary_fit$ok && primary_fit$beta["E"] > epsilon) {
  beta_full <- primary_fit$beta["E"]
  all_ids <- unique(panel$id)
  for (j in all_ids) {
    sub <- panel[panel$id != j,,drop=FALSE]
    cf <- coef_primary_only(sub)
    if (!cf$ok) {
      rel <- NA_real_
      reversal <- FALSE
      breach <- TRUE
    } else {
      rel <- abs(cf$beta - beta_full) / max(abs(beta_full),0.000001)
      reversal <- beta_full > epsilon && cf$beta < -epsilon
      breach <- reversal || rel >= 0.20
    }
    influence_rows[[length(influence_rows)+1]] <- data.frame(
      omitted_id=j,
      estimable=cf$ok,
      beta_minus_j=cf$beta,
      rank=cf$rank,
      rank_tol=cf$rank_tol,
      relative_change=rel,
      sign_reversal=reversal,
      breach=breach,
      reason=cf$reason,
      stringsAsFactors=FALSE
    )
  }
  influence <- do.call(rbind,influence_rows)
  n_breach <- sum(influence$breach)
  if (n_breach >= 5) {
    s5_status <- "FAIL"
  } else if (n_breach >= 1) {
    s5_status <- "PARTIAL_TRIGGER"
  } else {
    s5_status <- "PASS"
  }
} else {
  influence <- data.frame(
    omitted_id=character(),estimable=logical(),beta_minus_j=numeric(),
    rank=integer(),rank_tol=numeric(),relative_change=numeric(),
    sign_reversal=logical(),breach=logical(),reason=character()
  )
  n_breach <- NA_integer_
  s5_status <- if (primary_status == "FAIL") "NOT_RUN_PRIMARY_FAIL" else "FAIL"
}
write.csv(influence, file.path(out_dir, "S5_leave_one_organization_out.csv"), row.names=FALSE)

e0 <- as.numeric(panel$E_score == 0)
e2 <- as.numeric(panel$E_score == 2)
e3 <- as.numeric(panel$E_score == 3)
X6 <- cbind(
  Score0=demean_vec(e0,id),
  Score2=demean_vec(e2,id),
  Score3=demean_vec(e3,id),
  Wave2=demean_vec(w2,id),
  Wave3=demean_vec(w3,id),
  Wave4=demean_vec(w4,id)
)
s6_fit <- cluster_fit(
  X6,y3,id,
  c("Score0","Score2","Score3","Wave2","Wave3","Wave4"),
  "S6_nonlinear_exposure"
)
add_diag(s6_fit)

if (!s6_fit$ok) {
  s6_status <- "FAIL"
} else {
  neg_supported <- s6_fit$upper[c("Score2","Score3")] < -epsilon
  if (sum(neg_supported) == 2) {
    s6_status <- "FAIL"
  } else if (sum(neg_supported) == 1) {
    s6_status <- "PARTIAL_TRIGGER"
  } else {
    s6_status <- "PASS"
  }
}

s6_table <- if (s6_fit$ok) {
  data.frame(
    term=names(s6_fit$beta),
    estimate=as.numeric(s6_fit$beta),
    se=as.numeric(s6_fit$se),
    lower=as.numeric(s6_fit$lower),
    upper=as.numeric(s6_fit$upper)
  )
} else {
  data.frame(term="NON_ESTIMABLE",estimate=NA_real_,se=NA_real_,lower=NA_real_,upper=NA_real_)
}
write.csv(s6_table, file.path(out_dir, "S6_nonlinear_exposure.csv"), row.names=FALSE)

score_levels <- 0:4
ceiling_rows <- list()
for (wv in sort(unique(panel$wave))) {
  z <- panel[panel$wave == wv,,drop=FALSE]
  for (sc in score_levels) {
    n <- sum(z$R_score == sc)
    ceiling_rows[[length(ceiling_rows)+1]] <- data.frame(
      wave=wv,R_score=sc,count=n,share=n/nrow(z)
    )
  }
}
ceiling <- do.call(rbind,ceiling_rows)
write.csv(ceiling, file.path(out_dir, "ceiling_by_wave.csv"), row.names=FALSE)

transition_rows <- list()
for (w0 in 1:3) {
  q <- pairs[pairs$wave0 == w0 & pairs$wave1 == w0+1,,drop=FALSE]
  transition_rows[[length(transition_rows)+1]] <- data.frame(
    transition=paste0(w0,"_to_",w0+1),
    intervals=nrow(q),
    response_increase=sum(q$dR > epsilon),
    response_decrease=sum(q$dR < -epsilon),
    move_to_score4=sum(q$R0_score < 4 & q$R1_score == 4),
    move_away_from_score4=sum(q$R0_score == 4 & q$R1_score < 4)
  )
}
transition_diag <- do.call(rbind,transition_rows)
write.csv(transition_diag, file.path(out_dir, "adjacent_response_transition_diagnostics.csv"), row.names=FALSE)

sens_status <- c(
  S1=s1_status,S2=s2_status,S3=s3_status,
  S4=s4_status,S5=s5_status,S6=s6_status
)

final_status <- NA_character_
partial_triggers <- character()

if (primary_status == "FAIL") {
  final_status <- "LONGITUDINAL_FAIL"
} else {
  for (nm in names(sens_status)) {
    st <- sens_status[[nm]]
    if (identical(st,"FAIL")) {
      final_status <- "LONGITUDINAL_FAIL"
      break
    }
    if (identical(st,"PARTIAL_TRIGGER")) partial_triggers <- c(partial_triggers,nm)
  }
  if (is.na(final_status)) {
    if (primary_status == "PARTIAL") {
      final_status <- "LONGITUDINAL_PARTIAL"
    } else if (primary_status == "SUPPORT" && length(partial_triggers) == 0) {
      final_status <- "LONGITUDINAL_SUPPORT"
    } else {
      final_status <- "LONGITUDINAL_PARTIAL"
    }
  }
}

status_table <- data.frame(
  stage=c("Primary",names(sens_status)),
  status=c(primary_status,unname(sens_status)),
  stringsAsFactors=FALSE
)
write.csv(status_table, file.path(out_dir, "stage_statuses.csv"), row.names=FALSE)

diag_table <- do.call(rbind, diagnostic_rows)
write.csv(diag_table, file.path(out_dir, "model_rank_and_numerical_diagnostics.csv"), row.names=FALSE)

decision <- data.frame(
  final_status=final_status,
  primary_status=primary_status,
  primary_beta=if (primary_fit$ok) unname(primary_fit$beta["E"]) else NA_real_,
  primary_se=if (primary_fit$ok) unname(primary_fit$se["E"]) else NA_real_,
  primary_lower=if (primary_fit$ok) unname(primary_fit$lower["E"]) else NA_real_,
  primary_upper=if (primary_fit$ok) unname(primary_fit$upper["E"]) else NA_real_,
  partial_triggers=if (length(partial_triggers)) paste(partial_triggers,collapse=";") else "",
  S1=s1_status,S2=s2_status,S3=s3_status,S4=s4_status,S5=s5_status,S6=s6_status,
  stringsAsFactors=FALSE
)
write.csv(decision, file.path(out_dir, "longitudinal_decision.csv"), row.names=FALSE)

manifest <- c(
  paste("r_version",as.character(getRversion()),sep="="),
  paste("source_sha256",expected_source_sha256,sep="="),
  paste("raw_csv_sha256",expected_raw_csv_sha256,sep="="),
  paste("panel_rows",nrow(panel),sep="="),
  paste("organizations",length(unique(panel$id)),sep="="),
  paste("exposure_vary_orgs",n_e_var,sep="="),
  paste("response_vary_orgs",n_r_var,sep="="),
  paste("both_vary_orgs",n_both_var,sep="="),
  paste("adjacent_intervals",nrow(pairs),sep="="),
  paste("adjacent_exposure_vary_orgs",org_nonzero_de,sep="="),
  paste("adjacent_both_change_orgs",org_both_delta,sep="="),
  paste("R3_vary_orgs",r3_var_count,sep="="),
  paste("R3_vary_prop",format(r3_var_prop,digits=17),sep="="),
  paste("S5_breaches",n_breach,sep="="),
  paste("final_status",final_status,sep="=")
)
writeLines(manifest, file.path(out_dir, "execution_manifest.txt"))

cat("CSLS_LONGITUDINAL_EXECUTION_COMPLETE\n")
cat(paste0("final_status=",final_status,"\n"))
cat(paste0("primary_beta=",if(primary_fit$ok) format(primary_fit$beta["E"],digits=17) else "NA","\n"))
cat(paste0("primary_ci=[",if(primary_fit$ok) format(primary_fit$lower["E"],digits=17) else "NA",",",if(primary_fit$ok) format(primary_fit$upper["E"],digits=17) else "NA","]\n"))
cat(paste0("S1=",s1_status," S2=",s2_status," S3=",s3_status," S4=",s4_status," S5=",s5_status," S6=",s6_status,"\n"))

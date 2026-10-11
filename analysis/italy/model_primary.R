suppressPackageStartupMessages(library(survey))

if (as.character(getRversion()) != "4.6.1") stop("R version drift")
if (as.character(packageVersion("survey")) != "4.5") stop("survey version drift")

design <- readRDS("results/italy/build/italy_survey_design.rds")
d <- design$variables

if (nrow(d) != 16809) stop("analysis N drift")
if (degf(design) != 16808) stop("design degrees of freedom drift")

f0 <- response_breadth ~ footprint_breadth + size_factor + region_factor + macrosector_factor

m0 <- svyglm(
  f0,
  design=design,
  family=quasibinomial(link="logit"),
  rescale=FALSE
)

if (!isTRUE(m0$converged)) stop("IT_M0 did not converge")

coef_tab <- summary(m0)$coefficients
vc <- vcov(m0)
beta <- coef(m0)
df_resid <- m0$df.residual
w <- d$weight

make_x <- function(data_obj) {
  x <- model.matrix(delete.response(terms(f0)), data=data_obj)
  x[, names(beta), drop=FALSE]
}

weighted_prediction <- function(beta_value, data_obj) {
  x <- make_x(data_obj)
  p <- plogis(drop(x %*% beta_value))
  sum(w * p) / sum(w)
}

weighted_ame <- function(beta_value) {
  x <- make_x(d)
  p <- plogis(drop(x %*% beta_value))
  b <- unname(beta_value["footprint_breadth"])
  sum(w * b * p * (1-p)) / sum(w)
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
    g[j] <- (fn(bp)-fn(bm))/(2*h)
  }
  g
}

delta_result <- function(fn, beta_value, vcov_mat, df_resid) {
  est <- fn(beta_value)
  g <- num_grad(fn,beta_value)
  se <- sqrt(drop(t(g) %*% vcov_mat %*% g))
  crit <- qt(0.975,df=df_resid)
  c(estimate=est,se=se,lower=est-crit*se,upper=est+crit*se)
}

ame <- delta_result(weighted_ame,beta,vc,df_resid)

pred_rows <- list()
for (score in 0:4) {
  nd <- d
  nd$footprint_breadth <- score/4
  fn <- function(b) weighted_prediction(b,nd)
  dr <- delta_result(fn,beta,vc,df_resid)
  pred_rows[[length(pred_rows)+1]] <- data.frame(
    footprint_score=score,
    footprint_breadth=score/4,
    prediction=unname(dr["estimate"]),
    se=unname(dr["se"]),
    lower=unname(dr["lower"]),
    upper=unname(dr["upper"])
  )
}
pred_out <- do.call(rbind,pred_rows)

term <- "footprint_breadth"
est <- coef_tab[term,1]
se <- coef_tab[term,2]
pval <- coef_tab[term,4]
crit <- qt(0.975,df=df_resid)

primary <- data.frame(
  term=term,
  estimate=est,
  se=se,
  df=df_resid,
  p_value=pval,
  lower=est-crit*se,
  upper=est+crit*se
)

ame_out <- data.frame(
  estimand="weighted_average_marginal_effect",
  estimate=unname(ame["estimate"]),
  se=unname(ame["se"]),
  df=df_resid,
  lower=unname(ame["lower"]),
  upper=unname(ame["upper"])
)

fit <- fitted(m0)
checks <- c(
  n_match=nrow(d)==16809,
  design_degf_match=degf(design)==16808,
  converged=isTRUE(m0$converged),
  finite_coefficients=all(is.finite(coef_tab[,1])),
  finite_positive_se=all(is.finite(coef_tab[,2])) && all(coef_tab[,2] > 0),
  covariance_valid=all(is.finite(diag(vc))) && all(diag(vc) >= 0),
  fitted_unit_interval=all(is.finite(fit)) && all(fit >= 0 & fit <= 1),
  primary_finite=all(is.finite(unlist(primary[,c("estimate","se","p_value","lower","upper")]))),
  ame_finite=all(is.finite(unlist(ame_out[,c("estimate","se","lower","upper")]))),
  predictions_finite=all(is.finite(unlist(pred_out[,c("prediction","se","lower","upper")]))),
  predictions_unit_interval=all(pred_out$prediction >= 0 & pred_out$prediction <= 1)
)

write.csv(primary,"results/italy/models/italy_M0_primary_coefficient.csv",row.names=FALSE)
write.csv(ame_out,"results/italy/models/italy_M0_average_marginal_effect.csv",row.names=FALSE)
write.csv(pred_out,"results/italy/models/italy_M0_standardized_predictions.csv",row.names=FALSE)
write.csv(data.frame(term=rownames(coef_tab),estimate=coef_tab[,1],se=coef_tab[,2],statistic=coef_tab[,3],p_value=coef_tab[,4],row.names=NULL),
          "results/italy/models/italy_M0_full_coefficients.csv",row.names=FALSE)
writeLines(paste(names(checks),checks,sep="="),"results/italy/models/italy_M0_integrity_checks.txt")
writeLines(c(
  paste("n",nrow(d),sep="="),
  paste("design_degf",degf(design),sep="="),
  paste("model_df_residual",df_resid,sep="="),
  paste("converged",m0$converged,sep="="),
  paste("dispersion",summary(m0)$dispersion,sep="="),
  paste("r_version",as.character(getRversion()),sep="="),
  paste("survey_version",as.character(packageVersion("survey")),sep="=")
),"results/italy/models/italy_M0_model_manifest.txt")

saveRDS(m0,"results/italy/models/italy_M0_model.rds",version=3)

if (!all(checks)) quit(status=2)

cat("ITALY_M0_EXECUTION_PASS\n")
cat("estimate=",est,"\n",sep="")
cat("p=",pval,"\n",sep="")
cat("ame=",ame_out$estimate,"\n",sep="")

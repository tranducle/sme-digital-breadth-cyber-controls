suppressPackageStartupMessages(library(survey))

d <- readRDS("results/uk/build/uk_analysis_ready.rds")
if (nrow(d) != 2859) stop("S2 N drift")
if (!all(is.finite(d$weight)) || !all(d$weight > 0)) stop("raw weight invalid")

des <- svydesign(ids=~1, weights=~weight, data=d)

m0 <- svyglm(response_breadth ~ footprint_breadth + year_factor + size_medium + region_factor,
             design=des, family=quasibinomial(link="logit"), rescale=FALSE)
m1 <- svyglm(response_breadth ~ footprint_breadth * manage1 + year_factor + size_medium + region_factor,
             design=des, family=quasibinomial(link="logit"), rescale=FALSE)
m2 <- svyglm(response_breadth ~ footprint_breadth * manage2 + year_factor + size_medium + region_factor,
             design=des, family=quasibinomial(link="logit"), rescale=FALSE)

mods <- list(M0=m0,M1=m1,M2=m2)
if (!all(vapply(mods,function(x)isTRUE(x$converged),logical(1)))) stop("model convergence failure")

getrow <- function(model, term, label) {
  z <- summary(model)$coefficients
  data.frame(model=label, term=term, estimate=z[term,1], se=z[term,2], p_value=z[term,4])
}
out <- rbind(
  getrow(m0,"footprint_breadth","S2_M0"),
  getrow(m1,"footprint_breadth:manage1","S2_M1"),
  getrow(m2,"footprint_breadth:manage2","S2_M2")
)
out$holm_p <- NA_real_
idx <- out$model %in% c("S2_M1","S2_M2")
out$holm_p[idx] <- p.adjust(out$p_value[idx], method="holm")

checks <- c(
  n_match=nrow(d)==2859,
  all_converged=all(vapply(mods,function(x)isTRUE(x$converged),logical(1))),
  finite_estimates=all(is.finite(out$estimate)),
  finite_positive_se=all(is.finite(out$se) & out$se>0),
  finite_p=all(is.finite(out$p_value))
)

write.csv(out,"results/uk/sensitivity/S2_original_weights_key_results.csv",row.names=FALSE)
writeLines(paste(names(checks),checks,sep="="),"results/uk/sensitivity/S2_integrity_checks.txt")
writeLines(c(
 paste("n",nrow(d),sep="="),
 paste("weight_sum",sum(d$weight),sep="="),
 paste("weight_mean",mean(d$weight),sep="="),
 paste("r_version",as.character(getRversion()),sep="="),
 paste("survey_version",as.character(packageVersion("survey")),sep="=")
),"results/uk/sensitivity/S2_manifest.txt")

if (!all(checks)) quit(status=2)
cat("S2_EXECUTION_PASS\n")
print(out)

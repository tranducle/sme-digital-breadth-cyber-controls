suppressPackageStartupMessages(library(survey))

des <- readRDS("results/uk/build/uk_survey_design.rds")
d <- des$variables
if (nrow(d) != 2859) stop("S5 N drift")

m0 <- svyglm(response_score ~ footprint_breadth + year_factor + size_medium + region_factor,
             design=des, family=gaussian(), rescale=FALSE)
m1 <- svyglm(response_score ~ footprint_breadth * manage1 + year_factor + size_medium + region_factor,
             design=des, family=gaussian(), rescale=FALSE)
m2 <- svyglm(response_score ~ footprint_breadth * manage2 + year_factor + size_medium + region_factor,
             design=des, family=gaussian(), rescale=FALSE)

mods <- list(M0=m0,M1=m1,M2=m2)

getrow <- function(model,term,label) {
 z <- summary(model)$coefficients
 data.frame(model=label,term=term,estimate=z[term,1],se=z[term,2],p_value=z[term,4])
}

out <- rbind(
 getrow(m0,"footprint_breadth","S5_M0"),
 getrow(m1,"footprint_breadth:manage1","S5_M1"),
 getrow(m2,"footprint_breadth:manage2","S5_M2")
)
out$holm_p <- NA_real_
idx <- out$model %in% c("S5_M1","S5_M2")
out$holm_p[idx] <- p.adjust(out$p_value[idx],method="holm")

checks <- c(
 n_match=nrow(d)==2859,
 finite_estimates=all(is.finite(out$estimate)),
 finite_positive_se=all(is.finite(out$se)&out$se>0),
 finite_p=all(is.finite(out$p_value)),
 fitted_finite=all(vapply(mods,function(x)all(is.finite(fitted(x))),logical(1)))
)

write.csv(out,"results/uk/sensitivity/S5_linear_count_key_results.csv",row.names=FALSE)
writeLines(paste(names(checks),checks,sep="="),"results/uk/sensitivity/S5_integrity_checks.txt")

if (!all(checks)) quit(status=2)
cat("S5_EXECUTION_PASS\n")
print(out)

suppressPackageStartupMessages(library(survey))

des <- readRDS("results/uk/build/uk_survey_design.rds")
d <- des$variables
if (nrow(d) != 2859) stop("S6 N drift")

m <- svyglm(
 response_breadth ~ footprint_breadth * manage1 + footprint_breadth * manage2 +
 year_factor + size_medium + region_factor,
 design=des,
 family=quasibinomial(link="logit"),
 rescale=FALSE
)

if (!isTRUE(m$converged)) stop("S6 convergence failure")

z <- summary(m)$coefficients
terms <- c("footprint_breadth:manage1","footprint_breadth:manage2")
out <- data.frame(
 term=terms,
 estimate=z[terms,1],
 se=z[terms,2],
 p_value=z[terms,4],
 row.names=NULL
)
out$holm_p <- p.adjust(out$p_value,method="holm")

checks <- c(
 n_match=nrow(d)==2859,
 converged=isTRUE(m$converged),
 finite_estimates=all(is.finite(out$estimate)),
 finite_positive_se=all(is.finite(out$se)&out$se>0),
 finite_p=all(is.finite(out$p_value))
)

write.csv(out,"results/uk/sensitivity/S6_joint_organization_interactions.csv",row.names=FALSE)
write.csv(data.frame(term=rownames(z),estimate=z[,1],se=z[,2],statistic=z[,3],p_value=z[,4],row.names=NULL),
          "results/uk/sensitivity/S6_joint_full_coefficients.csv",row.names=FALSE)
writeLines(paste(names(checks),checks,sep="="),"results/uk/sensitivity/S6_integrity_checks.txt")

if (!all(checks)) quit(status=2)
cat("S6_EXECUTION_PASS\n")
print(out)

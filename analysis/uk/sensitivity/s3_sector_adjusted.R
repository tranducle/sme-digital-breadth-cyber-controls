suppressPackageStartupMessages(library(survey))

d <- readRDS("results/uk/build/uk_analysis_ready.rds")
x <- subset(d, year_calendar %in% 2019:2021)

if (nrow(x) != 1419) stop("S3 N drift")
if (any(is.na(x$sector_comb2))) stop("sector missing")

x$year_s3 <- factor(x$year_calendar, levels=2019:2021)
x$region_s3 <- factor(x$region_comb, levels=sort(unique(x$region_comb)))
x$sector_s3 <- factor(x$sector_comb2, levels=sort(unique(x$sector_comb2)))

des <- svydesign(ids=~1, weights=~w_norm, data=x)

base_controls <- "+ year_s3 + size_medium + region_s3 + sector_s3"

m0 <- svyglm(response_breadth ~ footprint_breadth + year_s3 + size_medium + region_s3 + sector_s3,
             design=des, family=quasibinomial(link="logit"), rescale=FALSE)
m1 <- svyglm(response_breadth ~ footprint_breadth * manage1 + year_s3 + size_medium + region_s3 + sector_s3,
             design=des, family=quasibinomial(link="logit"), rescale=FALSE)
m2 <- svyglm(response_breadth ~ footprint_breadth * manage2 + year_s3 + size_medium + region_s3 + sector_s3,
             design=des, family=quasibinomial(link="logit"), rescale=FALSE)

mods <- list(M0=m0,M1=m1,M2=m2)
if (!all(vapply(mods,function(z)isTRUE(z$converged),logical(1)))) stop("convergence failure")

getrow <- function(model,term,label) {
 z <- summary(model)$coefficients
 data.frame(model=label,term=term,estimate=z[term,1],se=z[term,2],p_value=z[term,4])
}
out <- rbind(
 getrow(m0,"footprint_breadth","S3_M0"),
 getrow(m1,"footprint_breadth:manage1","S3_M1"),
 getrow(m2,"footprint_breadth:manage2","S3_M2")
)
out$holm_p <- NA_real_
idx <- out$model %in% c("S3_M1","S3_M2")
out$holm_p[idx] <- p.adjust(out$p_value[idx], method="holm")

checks <- c(
 n_match=nrow(x)==1419,
 sector_complete=!any(is.na(x$sector_comb2)),
 all_converged=all(vapply(mods,function(z)isTRUE(z$converged),logical(1))),
 finite_estimates=all(is.finite(out$estimate)),
 finite_positive_se=all(is.finite(out$se)&out$se>0),
 finite_p=all(is.finite(out$p_value))
)

write.csv(out,"results/uk/sensitivity/S3_sector_adjusted_key_results.csv",row.names=FALSE)
writeLines(paste(names(checks),checks,sep="="),"results/uk/sensitivity/S3_integrity_checks.txt")
writeLines(c(
 paste("n",nrow(x),sep="="),
 paste("year_reference",levels(x$year_s3)[1],sep="="),
 paste("region_reference",levels(x$region_s3)[1],sep="="),
 paste("sector_reference",levels(x$sector_s3)[1],sep="="),
 paste("sector_levels",paste(levels(x$sector_s3),collapse=","),sep="=")
),"results/uk/sensitivity/S3_manifest.txt")

if (!all(checks)) quit(status=2)
cat("S3_EXECUTION_PASS\n")
print(out)

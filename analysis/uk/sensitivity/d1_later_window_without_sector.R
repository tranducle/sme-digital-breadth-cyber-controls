suppressPackageStartupMessages(library(survey))

d0 <- readRDS("results/uk/build/uk_analysis_ready.rds")
d <- subset(d0, year_calendar %in% 2019:2021)
if (nrow(d) != 1419) stop("D1 N drift")

d$year_d <- factor(d$year_calendar, levels=2019:2021)
d$region_d <- factor(d$region_comb, levels=sort(unique(d$region_comb)))

des <- svydesign(ids=~1, weights=~w_norm, data=d)

m0 <- svyglm(response_breadth ~ footprint_breadth + year_d + size_medium + region_d,
             design=des, family=quasibinomial(link="logit"), rescale=FALSE)
m1 <- svyglm(response_breadth ~ footprint_breadth * manage1 + year_d + size_medium + region_d,
             design=des, family=quasibinomial(link="logit"), rescale=FALSE)
m2 <- svyglm(response_breadth ~ footprint_breadth * manage2 + year_d + size_medium + region_d,
             design=des, family=quasibinomial(link="logit"), rescale=FALSE)

getrow <- function(model,term,label) {
 z <- summary(model)$coefficients
 data.frame(model=label,term=term,estimate=z[term,1],se=z[term,2],p_value=z[term,4])
}

out <- rbind(
 getrow(m0,"footprint_breadth","D1_M0"),
 getrow(m1,"footprint_breadth:manage1","D1_M1"),
 getrow(m2,"footprint_breadth:manage2","D1_M2")
)
out$holm_p <- NA_real_
out$holm_p[2:3] <- p.adjust(out$p_value[2:3],method="holm")

checks <- c(
 n_match=nrow(d)==1419,
 all_converged=all(c(m0$converged,m1$converged,m2$converged)),
 finite=all(is.finite(out$estimate)) && all(is.finite(out$se)) && all(is.finite(out$p_value))
)

write.csv(out,"results/uk/decomposition/D1_later_window_no_sector.csv",row.names=FALSE)
writeLines(paste(names(checks),checks,sep="="),"results/uk/decomposition/D1_checks.txt")

if (!all(checks)) quit(status=2)
cat("D1_EXECUTION_PASS\n")
print(out)

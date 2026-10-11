suppressPackageStartupMessages(library(survey))

d <- read.csv("results/uk/sensitivity/S1_unweighted_input.csv", stringsAsFactors=FALSE)
if (nrow(d) != 2860) stop("S1 N drift")

d$year_factor <- factor(d$year_calendar, levels=2017:2021)
d$region_factor <- factor(d$region_comb, levels=sort(unique(d$region_comb)))

des <- svydesign(ids=~1, weights=~one_weight, data=d)

m0 <- svyglm(response_breadth ~ footprint_breadth + year_factor + size_medium + region_factor,
             design=des, family=quasibinomial(link="logit"), rescale=FALSE)
m1 <- svyglm(response_breadth ~ footprint_breadth * manage1 + year_factor + size_medium + region_factor,
             design=des, family=quasibinomial(link="logit"), rescale=FALSE)
m2 <- svyglm(response_breadth ~ footprint_breadth * manage2 + year_factor + size_medium + region_factor,
             design=des, family=quasibinomial(link="logit"), rescale=FALSE)

mods <- list(M0=m0,M1=m1,M2=m2)
for (nm in names(mods)) {
  if (!isTRUE(mods[[nm]]$converged)) stop(paste(nm,"did not converge"))
}

getrow <- function(model, term, label) {
  z <- summary(model)$coefficients
  data.frame(model=label, term=term, estimate=z[term,1], se=z[term,2],
             p_value=z[term,4], stringsAsFactors=FALSE)
}

out <- rbind(
  getrow(m0,"footprint_breadth","S1_M0"),
  getrow(m1,"footprint_breadth:manage1","S1_M1"),
  getrow(m2,"footprint_breadth:manage2","S1_M2")
)

adj <- p.adjust(out$p_value[out$model %in% c("S1_M1","S1_M2")], method="holm")
out$holm_p <- NA_real_
out$holm_p[out$model %in% c("S1_M1","S1_M2")] <- adj

checks <- c(
  n_match=nrow(d)==2860,
  design_rows=nrow(des$variables)==2860,
  all_converged=all(vapply(mods,function(x)isTRUE(x$converged),logical(1))),
  finite_estimates=all(is.finite(out$estimate)),
  finite_se=all(is.finite(out$se) & out$se>0),
  finite_p=all(is.finite(out$p_value))
)

write.csv(out, "results/uk/sensitivity/S1_unweighted_key_results.csv", row.names=FALSE)
writeLines(paste(names(checks),checks,sep="="), "results/uk/sensitivity/S1_integrity_checks.txt")
writeLines(c(
  paste("n",nrow(d),sep="="),
  paste("design_degf",degf(des),sep="="),
  paste("r_version",as.character(getRversion()),sep="="),
  paste("survey_version",as.character(packageVersion("survey")),sep="=")
), "results/uk/sensitivity/S1_manifest.txt")

if (!all(checks)) quit(status=2)

cat("S1_EXECUTION_PASS\n")
print(out)

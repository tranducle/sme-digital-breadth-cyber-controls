suppressPackageStartupMessages(library(survey))

d0 <- readRDS("results/uk/build/uk_analysis_ready.rds")
if (nrow(d0) != 2859) stop("S7 N drift")

rows <- list()

for (omit in 2017:2021) {
  d <- subset(d0, year_calendar != omit)
  d$year_loo <- factor(d$year_calendar, levels=sort(unique(d$year_calendar)))
  d$region_loo <- factor(d$region_comb, levels=sort(unique(d$region_comb)))
  des <- svydesign(ids=~1, weights=~w_norm, data=d)

  m0 <- svyglm(response_breadth ~ footprint_breadth + year_loo + size_medium + region_loo,
               design=des, family=quasibinomial(link="logit"), rescale=FALSE)
  m1 <- svyglm(response_breadth ~ footprint_breadth * manage1 + year_loo + size_medium + region_loo,
               design=des, family=quasibinomial(link="logit"), rescale=FALSE)
  m2 <- svyglm(response_breadth ~ footprint_breadth * manage2 + year_loo + size_medium + region_loo,
               design=des, family=quasibinomial(link="logit"), rescale=FALSE)

  z0 <- summary(m0)$coefficients
  z1 <- summary(m1)$coefficients
  z2 <- summary(m2)$coefficients

  pmods <- c(z1["footprint_breadth:manage1",4], z2["footprint_breadth:manage2",4])
  padj <- p.adjust(pmods, method="holm")

  rows[[length(rows)+1]] <- data.frame(
    omitted_year=omit,
    n=nrow(d),
    m0_estimate=z0["footprint_breadth",1],
    m0_p=z0["footprint_breadth",4],
    m1_interaction=z1["footprint_breadth:manage1",1],
    m1_raw_p=pmods[1],
    m1_holm_p=padj[1],
    m2_interaction=z2["footprint_breadth:manage2",1],
    m2_raw_p=pmods[2],
    m2_holm_p=padj[2]
  )
}

out <- do.call(rbind,rows)

checks <- c(
 five_runs=nrow(out)==5,
 all_n_positive=all(out$n>0),
 finite_m0=all(is.finite(out$m0_estimate)) && all(is.finite(out$m0_p)),
 finite_m1=all(is.finite(out$m1_interaction)) && all(is.finite(out$m1_raw_p)),
 finite_m2=all(is.finite(out$m2_interaction)) && all(is.finite(out$m2_raw_p))
)

write.csv(out,"results/uk/sensitivity/S7_leave_one_wave_out.csv",row.names=FALSE)
writeLines(paste(names(checks),checks,sep="="),"results/uk/sensitivity/S7_integrity_checks.txt")

if (!all(checks)) quit(status=2)
cat("S7_EXECUTION_PASS\n")
print(out)

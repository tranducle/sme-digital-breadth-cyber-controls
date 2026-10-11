suppressPackageStartupMessages(library(survey))
des <- readRDS("results/italy/build/italy_survey_design.rds")
m <- svyglm(response_breadth ~ footprint_breadth + size_factor + region_factor + macrosector_factor + ict_domain_factor,
            design=des,family=quasibinomial(link="logit"),rescale=FALSE)
if (!isTRUE(m$converged)) stop("S6 convergence failure")
z <- summary(m)$coefficients
term <- "footprint_breadth"
out <- data.frame(model="IT_S6",term=term,estimate=z[term,1],se=z[term,2],p_value=z[term,4])
write.csv(out,"results/italy/models/italy_S6_ict_sector_adjusted.csv",row.names=FALSE)
cat("IT_S6_PASS\n")
print(out)

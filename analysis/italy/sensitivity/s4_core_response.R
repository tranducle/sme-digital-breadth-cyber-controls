suppressPackageStartupMessages(library(survey))
des <- readRDS("results/italy/build/italy_survey_design.rds")
d <- des$variables
if (!all(is.finite(d$security_core_breadth_6))) stop("S4 incomplete core response")
m <- svyglm(security_core_breadth_6 ~ footprint_breadth + size_factor + region_factor + macrosector_factor,
            design=des,family=quasibinomial(link="logit"),rescale=FALSE)
if (!isTRUE(m$converged)) stop("S4 convergence failure")
z <- summary(m)$coefficients
term <- "footprint_breadth"
out <- data.frame(model="IT_S4",term=term,estimate=z[term,1],se=z[term,2],p_value=z[term,4])
write.csv(out,"results/italy/models/italy_S4_core_response.csv",row.names=FALSE)
cat("IT_S4_PASS\n")
print(out)

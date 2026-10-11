suppressPackageStartupMessages(library(survey))
des <- readRDS("results/italy/build/italy_survey_design.rds")
d <- des$variables
if (!all(is.finite(d$digital_system_breadth_6))) stop("S3 incomplete broader Technology")
m <- svyglm(response_breadth ~ digital_system_breadth_6 + size_factor + region_factor + macrosector_factor,
            design=des,family=quasibinomial(link="logit"),rescale=FALSE)
if (!isTRUE(m$converged)) stop("S3 convergence failure")
z <- summary(m)$coefficients
term <- "digital_system_breadth_6"
out <- data.frame(model="IT_S3",term=term,estimate=z[term,1],se=z[term,2],p_value=z[term,4])
write.csv(out,"results/italy/models/italy_S3_broader_technology.csv",row.names=FALSE)
cat("IT_S3_PASS\n")
print(out)

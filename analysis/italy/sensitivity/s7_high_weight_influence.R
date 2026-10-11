suppressPackageStartupMessages(library(survey))
d <- readRDS("results/italy/build/italy_analysis_ready.rds")
thr <- 124.81329488119933
x <- subset(d, weight <= thr)
if (nrow(x) != 16640) stop(paste("S7 N drift",nrow(x)))
des <- svydesign(ids=~1,weights=~weight,data=x)
m <- svyglm(response_breadth ~ footprint_breadth + size_factor + region_factor + macrosector_factor,
            design=des,family=quasibinomial(link="logit"),rescale=FALSE)
if (!isTRUE(m$converged)) stop("S7 convergence failure")
z <- summary(m)$coefficients
term <- "footprint_breadth"
out <- data.frame(model="IT_S7",n=nrow(x),term=term,estimate=z[term,1],se=z[term,2],p_value=z[term,4])
write.csv(out,"results/italy/models/italy_S7_high_weight_influence.csv",row.names=FALSE)
cat("IT_S7_PASS\n")
print(out)

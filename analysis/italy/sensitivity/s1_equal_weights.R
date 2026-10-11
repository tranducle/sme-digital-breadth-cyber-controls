suppressPackageStartupMessages(library(survey))
d <- readRDS("results/italy/build/italy_analysis_ready.rds")
d$one_weight <- 1
des <- svydesign(ids=~1,weights=~one_weight,data=d)
m <- svyglm(response_breadth ~ footprint_breadth + size_factor + region_factor + macrosector_factor,
            design=des,family=quasibinomial(link="logit"),rescale=FALSE)
if (!isTRUE(m$converged)) stop("S1 convergence failure")
z <- summary(m)$coefficients
term <- "footprint_breadth"
out <- data.frame(
 model="IT_S1",
 term=term,
 estimate=z[term,1],
 se=z[term,2],
 p_value=z[term,4]
)
write.csv(out,"results/italy/models/italy_S1_equal_weight.csv",row.names=FALSE)
writeLines(c(
 paste("n",nrow(d),sep="="),
 paste("converged",m$converged,sep="="),
 paste("estimate",out$estimate,sep="="),
 paste("p",out$p_value,sep="=")
),"results/italy/models/italy_S1_manifest.txt")
cat("IT_S1_PASS\n")
print(out)

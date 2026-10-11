suppressPackageStartupMessages(library(survey))
d <- read.csv("results/italy/models/italy_S8_all_sizes_input.csv",stringsAsFactors=FALSE)
if (nrow(d)!=19804) stop("S8 N drift")
d$size3_factor <- factor(d$clad3,levels=c("cl1","cl2","cl3"))
d$region_factor <- factor(d$rip,levels=sort(unique(d$rip)))
d$macrosector_factor <- factor(d$mac,levels=sort(unique(d$mac)))
des <- svydesign(ids=~1,weights=~weight,data=d)
m <- svyglm(response_breadth ~ footprint_breadth + size3_factor + region_factor + macrosector_factor,
            design=des,family=quasibinomial(link="logit"),rescale=FALSE)
if (!isTRUE(m$converged)) stop("S8 convergence failure")
z <- summary(m)$coefficients
term <- "footprint_breadth"
out <- data.frame(model="IT_S8",n=nrow(d),term=term,estimate=z[term,1],se=z[term,2],p_value=z[term,4])
write.csv(out,"results/italy/models/italy_S8_large_enterprise_scope.csv",row.names=FALSE)
cat("IT_S8_PASS\n")
print(out)

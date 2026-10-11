suppressPackageStartupMessages(library(survey))
des <- readRDS("results/italy/build/italy_survey_design.rds")
m <- svyglm(response_score ~ footprint_breadth + size_factor + region_factor + macrosector_factor,
            design=des,family=gaussian(),rescale=FALSE)
z <- summary(m)$coefficients
term <- "footprint_breadth"
out <- data.frame(model="IT_S5",term=term,estimate=z[term,1],se=z[term,2],p_value=z[term,4])
write.csv(out,"results/italy/models/italy_S5_linear_response.csv",row.names=FALSE)
cat("IT_S5_PASS\n")
print(out)

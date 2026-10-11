suppressPackageStartupMessages(library(survey))
des <- readRDS("results/italy/build/italy_survey_design.rds")
des <- update(des, footprint_factor=factor(footprint_score,levels=0:4))
d <- des$variables

m <- svyglm(response_breadth ~ footprint_factor + size_factor + region_factor + macrosector_factor,
            design=des,family=quasibinomial(link="logit"),rescale=FALSE)
if (!isTRUE(m$converged)) stop("S2 convergence failure")

tt <- regTermTest(m,~footprint_factor)

beta <- coef(m)
w <- d$weight
terms_no_y <- delete.response(terms(m))
pred_rows <- list()
for (score in 0:4) {
  nd <- d
  nd$footprint_factor <- factor(score,levels=0:4)
  X <- model.matrix(terms_no_y,data=nd)
  X <- X[,names(beta),drop=FALSE]
  pr <- plogis(drop(X %*% beta))
  pred_rows[[length(pred_rows)+1]] <- data.frame(
    footprint_score=score,
    standardized_prediction=sum(w*pr)/sum(w)
  )
}
pred <- do.call(rbind,pred_rows)

out <- data.frame(
  model="IT_S2",
  term="footprint_factor",
  numerator_df=tt$df,
  denominator_df=tt$ddf,
  F_stat=as.numeric(tt$Ftest),
  p_value=as.numeric(tt$p)
)

write.csv(out,"results/italy/models/italy_S2_categorical_omnibus.csv",row.names=FALSE)
write.csv(pred,"results/italy/models/italy_S2_categorical_predictions.csv",row.names=FALSE)
cat("IT_S2_PASS\n")
print(out)
print(pred)

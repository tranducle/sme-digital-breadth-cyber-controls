suppressPackageStartupMessages(library(survey))

des <- readRDS("results/uk/build/uk_survey_design.rds")
d <- des$variables
if (nrow(d) != 2859) stop("S4 N drift")

des <- update(des, footprint_factor=factor(footprint_score, levels=0:4))
d <- des$variables

m0 <- svyglm(response_breadth ~ footprint_factor + year_factor + size_medium + region_factor,
             design=des, family=quasibinomial(link="logit"), rescale=FALSE)
m1 <- svyglm(response_breadth ~ footprint_factor * manage1 + year_factor + size_medium + region_factor,
             design=des, family=quasibinomial(link="logit"), rescale=FALSE)
m2 <- svyglm(response_breadth ~ footprint_factor * manage2 + year_factor + size_medium + region_factor,
             design=des, family=quasibinomial(link="logit"), rescale=FALSE)

if (!all(c(m0$converged,m1$converged,m2$converged))) stop("S4 convergence failure")

t0 <- regTermTest(m0, ~footprint_factor)
t1 <- regTermTest(m1, ~footprint_factor:manage1)
t2 <- regTermTest(m2, ~footprint_factor:manage2)

tests <- data.frame(
 model=c("S4_M0","S4_M1","S4_M2"),
 term=c("footprint_factor","footprint_factor:manage1","footprint_factor:manage2"),
 numerator_df=c(t0$df,t1$df,t2$df),
 denominator_df=c(t0$ddf,t1$ddf,t2$ddf),
 F_stat=c(as.numeric(t0$Ftest),as.numeric(t1$Ftest),as.numeric(t2$Ftest)),
 p_value=c(as.numeric(t0$p),as.numeric(t1$p),as.numeric(t2$p))
)
tests$holm_p <- NA_real_
tests$holm_p[2:3] <- p.adjust(tests$p_value[2:3], method="holm")

weighted_standardized <- function(model, moderator_name=NULL, moderator_value=NULL) {
  beta <- coef(model)
  terms_no_y <- delete.response(terms(model))
  w <- d$w_norm
  rows <- list()
  for (score in 0:4) {
    nd <- d
    nd$footprint_factor <- factor(score, levels=0:4)
    if (!is.null(moderator_name)) nd[[moderator_name]] <- moderator_value
    X <- model.matrix(terms_no_y, data=nd)
    X <- X[,names(beta),drop=FALSE]
    pred <- plogis(drop(X %*% beta))
    rows[[length(rows)+1]] <- data.frame(
      footprint_score=score,
      prediction=sum(w*pred)/sum(w)
    )
  }
  do.call(rbind,rows)
}

pred0 <- weighted_standardized(m0)
pred0$model <- "S4_M0"
pred0$moderator <- NA_character_
pred0$state <- NA_integer_

pred1_0 <- weighted_standardized(m1,"manage1",0)
pred1_0$model <- "S4_M1"
pred1_0$moderator <- "manage1"
pred1_0$state <- 0L
pred1_1 <- weighted_standardized(m1,"manage1",1)
pred1_1$model <- "S4_M1"
pred1_1$moderator <- "manage1"
pred1_1$state <- 1L

pred2_0 <- weighted_standardized(m2,"manage2",0)
pred2_0$model <- "S4_M2"
pred2_0$moderator <- "manage2"
pred2_0$state <- 0L
pred2_1 <- weighted_standardized(m2,"manage2",1)
pred2_1$model <- "S4_M2"
pred2_1$moderator <- "manage2"
pred2_1$state <- 1L

pred <- rbind(pred0,pred1_0,pred1_1,pred2_0,pred2_1)

checks <- c(
 n_match=nrow(d)==2859,
 all_converged=all(c(m0$converged,m1$converged,m2$converged)),
 all_tests_finite=all(is.finite(tests$F_stat)) && all(is.finite(tests$p_value)),
 predictions_finite=all(is.finite(pred$prediction)),
 predictions_unit_interval=all(pred$prediction>=0 & pred$prediction<=1)
)

write.csv(tests,"results/uk/sensitivity/S4_categorical_omnibus_tests.csv",row.names=FALSE)
write.csv(pred,"results/uk/sensitivity/S4_categorical_standardized_predictions.csv",row.names=FALSE)
writeLines(paste(names(checks),checks,sep="="),"results/uk/sensitivity/S4_integrity_checks.txt")

if (!all(checks)) quit(status=2)
cat("S4_EXECUTION_PASS\n")
print(tests)
cat("M1 predictions\n")
print(subset(pred,model=="S4_M1"))

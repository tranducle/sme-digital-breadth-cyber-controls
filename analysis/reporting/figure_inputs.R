# Assemble the plotted values for the manuscript figures.
# Every value except the digital-breadth score distribution is copied from an
# analysis output; the distribution is computed here from the survey designs. Writes one CSV per figure to results/figure_inputs.
# Run from the repository root.

suppressPackageStartupMessages(library(survey))

out_dir <- "results/figure_inputs"
exp_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
rd <- function(...) read.csv(file.path(exp_dir, ...))

# Figure 2: predictions by digital-breadth score, plus score distribution.
uk_des <- readRDS(file.path(exp_dir,
  "uk/build/uk_survey_design.rds"))
it_des <- readRDS(file.path(exp_dir,
  "italy/build/italy_survey_design.rds"))
if (nrow(uk_des$variables) != 2859) stop("UK N drift")
if (nrow(it_des$variables) != 16809) stop("Italy N drift")

score_share <- function(des) {
  m <- svymean(~factor(footprint_score), des)
  lv <- as.integer(sub("^factor\\(footprint_score\\)", "", names(coef(m))))
  data.frame(footprint_score = lv,
             n = as.vector(table(des$variables$footprint_score)[as.character(lv)]),
             weighted_share = as.vector(coef(m)))
}

uk_m0 <- rd("uk/primary/M0_standardized_predictions.csv")
it_m0 <- rd("italy/models/italy_M0_standardized_predictions.csv")
uk_cat <- subset(rd("uk/sensitivity/S4_categorical_standardized_predictions.csv"),
                 model == "S4_M0")
it_cat <- rd("italy/models/italy_S2_categorical_predictions.csv")

fig2 <- rbind(
  data.frame(sample = "UK", uk_m0[, c("footprint_score", "prediction", "lower", "upper")],
             categorical = uk_cat$prediction[match(uk_m0$footprint_score, uk_cat$footprint_score)],
             score_share(uk_des)[, c("n", "weighted_share")]),
  data.frame(sample = "Italy", it_m0[, c("footprint_score", "prediction", "lower", "upper")],
             categorical = it_cat$standardized_prediction[match(it_m0$footprint_score, it_cat$footprint_score)],
             score_share(it_des)[, c("n", "weighted_share")])
)
if (any(abs(tapply(fig2$weighted_share, fig2$sample, sum) - 1) > 1e-8)) stop("Share does not sum to 1")
write.csv(fig2, file.path(out_dir, "figure2_predictions_by_score.csv"), row.names = FALSE)

# Figure 3: CSLS response distribution by wave and RQ3 estimates.
csls <- "panel"
write.csv(rd(csls, "ceiling_by_wave.csv"), file.path(out_dir, "figure3a_response_by_wave.csv"), row.names = FALSE)
pick <- function(file, term, label) {
  x <- rd(csls, file); x <- x[x$term == term, ]
  data.frame(specification = label, estimate = x$estimate, lower = x$lower, upper = x$upper)
}
loo <- rd(csls, "S5_leave_one_organization_out.csv")
loo <- loo[loo$estimable, ]
fig3b <- rbind(
  pick("primary_coefficients.csv", "E", "Primary fixed-effects model"),
  pick("S1_first_difference.csv", "dE", "S1 First differences"),
  pick("S2_monitoring_excluded.csv", "E", "S2 Monitoring excluded"),
  pick("S4_mode_adjustment.csv", "E", "S4 Survey-mode adjustment"),
  data.frame(specification = "S5 Leave one organization out",
             estimate = NA, lower = min(loo$beta_minus_j), upper = max(loo$beta_minus_j))
)
write.csv(fig3b, file.path(out_dir, "figure3b_rq3_estimates.csv"), row.names = FALSE)

# Figure 4: next-wave contrasts, primary and the two timing checks.
st_contrasts <- rd("next_wave/state_transition/state_transition_contrasts.csv")
cfg_contrasts <- rd("next_wave/configuration/configuration_contrasts.csv")
spec_st <- c(PRIMARY = "All adjacent pairs", S4 = "Wave 3 to 4 only", S5 = "Two-wave gap")
spec_cfg <- c(primary = "All adjacent pairs", S4_strict_wave3_to4 = "Wave 3 to 4 only",
             S5_two_wave_gap = "Two-wave gap")
fig4 <- do.call(rbind, c(
  lapply(c("O", "S"), function(q) {
    x <- st_contrasts[st_contrasts$contrast == q & st_contrasts$spec %in% names(spec_st), ]
    data.frame(quantity = q, specification = spec_st[x$spec],
               estimate = x$estimate, lower = x$lcl, upper = x$ucl)
  }),
  list({
    x <- cfg_contrasts[cfg_contrasts$protection_state == "consequence" & cfg_contrasts$spec %in% names(spec_cfg), ]
    data.frame(quantity = "C", specification = spec_cfg[x$spec],
               estimate = x$estimate, lower = x$lower, upper = x$upper)
  })
))
write.csv(fig4, file.path(out_dir, "figure4_next_wave_contrasts.csv"), row.names = FALSE)

print(fig2); print(fig3b); print(fig4)
cat("FIGURE_INPUTS_PASS\n")

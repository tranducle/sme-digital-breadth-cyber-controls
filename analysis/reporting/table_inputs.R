# Compute descriptive inputs for manuscript tables (sample composition,
# weighted item prevalence, mean breadth) and refit UK leave-one-wave-out
# M0 models to obtain standard errors that S7 did not store.
# Reads only analysis-ready artifacts; writes CSVs to results/table_inputs.
# Run from the repository root.

suppressPackageStartupMessages(library(survey))

out_dir <- "results/table_inputs"
paper_root <- "."
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

uk_des <- readRDS(file.path(paper_root,
  "results/uk/build/uk_survey_design.rds"))
it_des <- readRDS(file.path(paper_root,
  "results/italy/build/italy_survey_design.rds"))
csls <- read.csv(file.path(paper_root,
  "results/panel/analysis_ready_panel.csv"))

if (nrow(uk_des$variables) != 2859) stop("UK N drift")
if (nrow(it_des$variables) != 16809) stop("Italy N drift")
if (nrow(csls) != 783 || length(unique(csls$id)) != 318) stop("CSLS N drift")

wshare <- function(des, var) {
  f <- as.formula(paste0("~factor(", var, ")"))
  m <- svymean(f, des)
  lv <- sub(paste0("^factor\\(", var, "\\)"), "", names(coef(m)))
  data.frame(variable = var, level = lv,
             n = as.vector(table(des$variables[[var]])[lv]),
             weighted_share = as.vector(coef(m)))
}

# 1. Sample composition
comp <- rbind(
  cbind(sample = "UK", wshare(uk_des, "year_calendar")),
  cbind(sample = "UK", wshare(uk_des, "sizeb")),
  cbind(sample = "UK", wshare(uk_des, "region_comb")),
  cbind(sample = "Italy", wshare(it_des, "clad3")),
  cbind(sample = "Italy", wshare(it_des, "rip")),
  cbind(sample = "Italy", wshare(it_des, "mac"))
)
csls_wave <- aggregate(id ~ wave, data = csls, FUN = length)
comp <- rbind(comp, data.frame(sample = "CSLS", variable = "wave",
  level = csls_wave$wave, n = csls_wave$id,
  weighted_share = csls_wave$id / nrow(csls)))
write.csv(comp, file.path(out_dir, "table_sample_composition.csv"), row.names = FALSE)

# 2. Item prevalence (weighted for national samples, unweighted for CSLS)
uk_bin <- function(x) as.numeric(x == 1)
uk_des <- update(uk_des,
  i_online4 = uk_bin(online4), i_online5 = uk_bin(online5),
  i_online9 = uk_bin(online9), i_mobile = as.numeric(mobile_bin == 1),
  i_rules2 = uk_bin(rules2), i_rules4 = uk_bin(rules4),
  i_rules5 = uk_bin(rules5), i_rules7 = uk_bin(rules7))
uk_items <- c("i_online4", "i_online5", "i_online9", "i_mobile",
              "i_rules2", "i_rules4", "i_rules5", "i_rules7")
# Consistency check: item indicators must reproduce the stored scores.
v <- uk_des$variables
if (any(rowSums(v[, uk_items[1:4]]) != v$footprint_score)) stop("UK footprint item mismatch")
if (any(rowSums(v[, uk_items[5:8]]) != v$response_score)) stop("UK response item mismatch")

it_items <- c("c5_num", "c7_num", "C8a_routed", "C10a_num",
              paste0("e1", letters[1:10], "_num"))
w <- it_des$variables
if (any(rowSums(w[, it_items[1:4]]) != w$footprint_score)) stop("Italy footprint item mismatch")
if (any(rowSums(w[, it_items[5:14]]) != w$response_score)) stop("Italy response item mismatch")

prev_svy <- function(des, items, sample) {
  m <- svymean(as.formula(paste("~", paste(items, collapse = "+"))), des)
  data.frame(sample = sample, item = items,
             unweighted = colMeans(des$variables[, items]),
             weighted = as.vector(coef(m)), se = as.vector(SE(m)))
}
csls_items <- c("onlinea_bin", "onlineb_bin", "devices_bin",
                "rulesb_bin", "rulesf_bin", "rulesh_bin", "rulesi_bin")
prev <- rbind(
  prev_svy(uk_des, uk_items, "UK"),
  prev_svy(it_des, it_items, "Italy"),
  data.frame(sample = "CSLS", item = csls_items,
             unweighted = colMeans(csls[, csls_items]),
             weighted = NA, se = NA)
)
write.csv(prev, file.path(out_dir, "table_item_prevalence.csv"), row.names = FALSE)

# 3. Mean breadth and share at maximum response
brd <- function(des, sample, max_resp) {
  m <- svymean(~footprint_breadth + response_breadth, des)
  at_max <- svymean(~I(as.numeric(response_score == max_resp)), des)
  data.frame(sample = sample,
             mean_footprint_breadth = coef(m)[1], mean_response_breadth = coef(m)[2],
             share_response_at_max = coef(at_max)[1])
}
breadth <- rbind(brd(uk_des, "UK", 4), brd(it_des, "Italy", 10),
  data.frame(sample = "CSLS", mean_footprint_breadth = mean(csls$E),
             mean_response_breadth = mean(csls$R),
             share_response_at_max = mean(csls$R_score == 4)))
write.csv(breadth, file.path(out_dir, "table_mean_breadth.csv"), row.names = FALSE)

# 4. UK leave-one-wave-out refit with standard errors (same formula as S7)
d0 <- uk_des$variables
lowo <- do.call(rbind, lapply(2017:2021, function(omit) {
  d <- subset(d0, year_calendar != omit)
  d$year_loo <- factor(d$year_calendar, levels = sort(unique(d$year_calendar)))
  d$region_loo <- factor(d$region_comb, levels = sort(unique(d$region_comb)))
  des <- svydesign(ids = ~1, weights = ~w_norm, data = d)
  m0 <- svyglm(response_breadth ~ footprint_breadth + year_loo + size_medium + region_loo,
               design = des, family = quasibinomial(link = "logit"), rescale = FALSE)
  z <- summary(m0)$coefficients["footprint_breadth", ]
  dfr <- m0$df.residual
  data.frame(omitted_year = omit, n = nrow(d), estimate = z[1], se = z[2],
             df = dfr, lower = z[1] - qt(0.975, dfr) * z[2],
             upper = z[1] + qt(0.975, dfr) * z[2], p_value = z[4])
}))
write.csv(lowo, file.path(out_dir, "table_uk_lowo_with_se.csv"), row.names = FALSE)

print(comp); print(prev); print(breadth); print(lowo)
cat("MANUSCRIPT_TABLE_INPUTS_PASS\n")

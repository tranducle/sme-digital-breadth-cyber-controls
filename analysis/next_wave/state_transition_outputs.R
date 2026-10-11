# State-transition output writer. Every artifact is written whether the run completed or failed.

st_num <- function(x) if (length(x) == 1 && !is.null(x)) as.numeric(x) else NA_real_
st_or_none <- function(x) if (length(x)) paste(x, collapse = ",") else "none"

st_summary_row <- function(r) {
  g <- function(term, field) if (isTRUE(r$ok)) st_num(r$coefs[r$coefs$term == term, field]) else NA_real_
  k <- function(name, field) if (isTRUE(r$ok)) st_num(r$contrasts[r$contrasts$contrast == name, field]) else NA_real_
  st <- if (is.null(r$states)) c(NA_real_, NA_real_) else r$states
  data.frame(
    spec = r$key, label = r$label, estimable = isTRUE(r$ok),
    failure_reasons = paste(r$reasons, collapse = "; "),
    separation_warning = isTRUE(r$separation), n_warnings = length(r$warnings),
    N = st_num(r$n), G = st_num(r$G), n_incident_free = st_num(r$n_it0),
    n_incident_positive = st_num(r$n_it1), events = st_num(r$events),
    events_incident_free = st_num(r$events_it0),
    protection_low_state = st[1], protection_high_state = st[2],
    cr1_multiplier = if (isTRUE(r$ok)) r$vc$cr1 else NA_real_,
    t_df = if (isTRUE(r$ok)) r$vc$G - 1 else NA_real_,
    betaR = g("Rc", "estimate"), betaR_se = g("Rc", "se"), betaR_lcl = g("Rc", "lcl"), betaR_ucl = g("Rc", "ucl"),
    betaIR = g("It:Rc", "estimate"), betaIR_se = g("It:Rc", "se"),
    betaIR_lcl = g("It:Rc", "lcl"), betaIR_ucl = g("It:Rc", "ucl"),
    betaIE = g("It:Ec", "estimate"), betaIE_se = g("It:Ec", "se"),
    O = k("O", "estimate"), O_se = k("O", "se"), O_lcl = k("O", "lcl"), O_ucl = k("O", "ucl"),
    P = k("P", "estimate"), P_se = k("P", "se"), P_lcl = k("P", "lcl"), P_ucl = k("P", "ucl"),
    O_common = k("O_common", "estimate"), P_common = k("P_common", "estimate"),
    S = k("S", "estimate"), S_se = k("S", "se"), S_lcl = k("S", "lcl"), S_ucl = k("S", "ucl"),
    contradictory = isTRUE(r$ok) && st_is_contradictory(g("Rc", "lcl"), k("O", "estimate")),
    directional_onset = isTRUE(r$ok) && g("Rc", "estimate") < 0 && k("O", "estimate") > 0,
    descriptive_O_ci_below_0 = isTRUE(r$ok) && k("O", "ucl") < 0,
    stringsAsFactors = FALSE)
}

st_stack <- function(results, field) {
  parts <- lapply(results, function(r) if (isTRUE(r$ok)) cbind(spec = r$key, r[[field]]) else NULL)
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(data.frame(spec = character(0)))
  do.call(rbind, parts)
}

st_conditions_frame <- function(dec, results) {
  if (is.null(dec$conditions)) {
    return(data.frame(condition = "NOT_EVALUATED", pass = FALSE, evidence = dec$failure_reason))
  }
  p <- results$PRIMARY
  g <- function(term, field) signif(p$coefs[p$coefs$term == term, field], 8)
  k <- function(name, field) signif(p$contrasts[p$contrasts$contrast == name, field], 8)
  evidence <- c(
    sprintf("betaR=%s", g("Rc", "estimate")),
    sprintf("betaR_ucl=%s", g("Rc", "ucl")),
    sprintf("O=%s threshold=%s", k("O", "estimate"), ST_THRESHOLD_O),
    sprintf("O_lcl=%s", k("O", "lcl")),
    sprintf("betaIR_ci=[%s, %s]", g("It:Rc", "lcl"), g("It:Rc", "ucl")),
    sprintf("abs(S)=%s threshold=%s", abs(k("S", "estimate")), ST_THRESHOLD_S),
    sprintf("S_ci=[%s, %s]", k("S", "lcl"), k("S", "ucl")),
    sprintf("sign(betaIR)=%s sign(S)=%s", sign(g("It:Rc", "estimate")), sign(k("S", "estimate"))),
    sprintf("contradictory_specs=%s", st_or_none(dec$contradictory_specs))
  )
  data.frame(condition = names(dec$conditions), pass = unname(dec$conditions), evidence = evidence)
}

st_retention_frame <- function(ret) {
  if (!isTRUE(ret$ok)) return(data.frame(metric = "unavailable", probability = NA, normalized_weight = NA))
  q <- function(x) c(min(x), unname(quantile(x, 0.05)), median(x), unname(quantile(x, 0.95)), max(x), mean(x))
  data.frame(metric = c("min", "p05", "median", "p95", "max", "mean"),
             probability = q(ret$p), normalized_weight = q(ret$w))
}

st_decision_list <- function(state) {
  d <- state$decision
  list(
    model = "state_transition", mode = state$args$mode, run_status = state$status,
    verdict = d$verdict, integrity_ok = isTRUE(d$integrity_ok), failure_reason = d$failure_reason,
    support_conditions = if (is.null(d$conditions)) NULL else as.list(d$conditions),
    contradictory_specs = I(d$contradictory_specs), non_estimable_specs = I(d$non_estimable_specs),
    influence_triggered = isTRUE(d$influence_triggered), influence_class = d$influence_class,
    influence_breaches = d$influence_breaches,
    timing_onset_language_allowed = isTRUE(d$timing_onset_language_allowed),
    state_dependence_temporally_replicated = isTRUE(d$state_dependence_temporally_replicated),
    claim_status = d$claim_status,
    primary_separation_warning = isTRUE(state$results$PRIMARY$separation),
    sensitivity_separation_warnings = I(as.character(names(Filter(function(r) isTRUE(r$separation), state$results[-1])))),
    thresholds = list(O = ST_THRESHOLD_O, S = ST_THRESHOLD_S, sign_tol = ST_SIGN_TOL,
                      influence_rel = ST_INFLUENCE_REL, influence_fail_count = ST_INFLUENCE_FAIL_COUNT)
  )
}

st_write_outputs <- function(state, out) {
  wcsv <- function(x, name) write.csv(x, out(name), row.names = FALSE)
  res <- state$results
  checks <- if (is.null(state$checks)) data.frame(check = "pre_estimation_checks", observed = NA,
                                                 expected = NA, pass = FALSE) else state$checks
  wcsv(checks, "state_transition_integrity_checks.csv")
  wcsv(if (length(res)) do.call(rbind, lapply(res, st_summary_row)) else data.frame(spec = character(0)),
       "state_transition_specification_summary.csv")
  wcsv(st_stack(res, "coefs"), "state_transition_all_coefficients.csv")
  wcsv(if (isTRUE(res$PRIMARY$ok)) res$PRIMARY$coefs else data.frame(term = character(0)),
       "state_transition_primary_coefficients.csv")
  wcsv(st_stack(res, "contrasts"), "state_transition_contrasts.csv")
  wcsv(st_conditions_frame(state$decision, res), "state_transition_support_conditions.csv")
  infl <- if (is.null(state$influence)) data.frame(status = "NOT_TRIGGERED",
    reason = "influence diagnostic runs only when the primary verdict is SUPPORT or PARTIAL") else state$influence
  wcsv(infl, "state_transition_influence_diagnostic.csv")
  wcsv(st_retention_frame(state$retention), "state_transition_retention_weight_diagnostics.csv")
  warn <- do.call(rbind, lapply(res, function(r) if (length(r$warnings)) data.frame(spec = r$key, warning = r$warnings)))
  wcsv(if (is.null(warn)) data.frame(spec = character(0), warning = character(0)) else warn, "state_transition_warnings.csv")
  dl <- st_decision_list(state)
  jsonlite::write_json(dl, out("state_transition_decision.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
  objs <- lapply(res, function(r) list(fit = r$fit, vc = r$vc, sample = r$sample, ok = r$ok))
  saveRDS(list(results = objs, retention = state$retention, influence = state$influence), out("state_transition_model_objects.rds"))
  meta <- list(
    mode = state$args$mode, started = state$started, finished = state$finished, run_status = state$status,
    input_dir = state$args$input_dir, output_dir = state$args$output_dir,
    input_sha256 = as.list(state$hashes), script_sha256 = as.list(state$script_hashes),
    r_version = R.version.string,
    packages = list(digest = as.character(packageVersion("digest")),
                    jsonlite = as.character(packageVersion("jsonlite"))),
    centering = lapply(Filter(function(r) !is.null(r$sample), res),
                       function(r) list(r_center = r$sample$r_center, e_center = r$sample$e_center)),
    platform = R.version$platform
  )
  jsonlite::write_json(meta, out("state_transition_run_metadata.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  lines <- c(
    "State-transition decision record",
    paste("mode:", dl$mode), paste("run_status:", dl$run_status), paste("verdict:", dl$verdict),
    paste("integrity_ok:", dl$integrity_ok), paste("failure_reason:", dl$failure_reason),
    paste("contradictory_specs:", st_or_none(dl$contradictory_specs)),
    paste("non_estimable_specs:", st_or_none(dl$non_estimable_specs)),
    paste("influence_class:", dl$influence_class), paste("influence_breaches:", dl$influence_breaches),
    paste("timing_onset_language_allowed:", dl$timing_onset_language_allowed),
    paste("state_dependence_temporally_replicated:", dl$state_dependence_temporally_replicated),
    paste("claim_status:", dl$claim_status),
    paste("finished:", state$finished)
  )
  if (!is.null(dl$support_conditions)) {
    lines <- c(lines, paste0("condition ", names(dl$support_conditions), ": ", unlist(dl$support_conditions)))
  }
  writeLines(lines, out("state_transition_decision.txt"))
}

# Last-resort record when st_write_outputs itself fails: the decision files still state FAIL and
# the output set stays incomplete.
st_write_fallback <- function(state, out, err) {
  dl <- list(model = "state_transition", mode = state$args$mode, run_status = "OUTPUT_WRITE_FAILURE",
             verdict = ST_VERDICT_FAIL, integrity_ok = FALSE,
             failure_reason = paste0("output writer failure: ", err, " (run status ", state$status, ")"),
             claim_status = "NO_CLAIM")
  jsonlite::write_json(dl, out("state_transition_decision.json"), auto_unbox = TRUE, pretty = TRUE)
  writeLines(c("State-transition decision record", paste0(names(dl), ": ", unlist(dl))), out("state_transition_decision.txt"))
}

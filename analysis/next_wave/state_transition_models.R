# State-transition model driver.
# Usage: Rscript analysis/next_wave/state_transition_models.R --input-dir DIR --output-dir DIR
# The input files must match the SHA-256 values in ST_EXPECTED_INPUTS (the files written by
# build_inputs.py). Exit codes: 0 completed, 1 refused because an input hash differs (nothing
# written), 3 integrity or writer failure (the decision record is still written).
options(stringsAsFactors = FALSE, warn = 1)

st_script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
st_script_dir <- dirname(st_script_path)
source(file.path(st_script_dir, "state_transition_functions.R"))
source(file.path(st_script_dir, "state_transition_outputs.R"))

ST_EXPECTED_INPUTS <- c(
  configuration_adjacent_pairs.csv = "c67e973444186d248783f61c355b092e31c139df8f243ebb8eb2dd780060e725",
  configuration_two_wave_gap_pairs.csv = "3fc47646531167a90d82f9741549602d200cecfe553b3fe08bdb836534872f3d",
  configuration_retention_eligibility.csv = "3ae45c4dfdd2dfefbf28aac585c9a0e9625720a26f887db52c8504e1ae56c1a0"
)
ST_SPECS <- list(
  PRIMARY = list(label = "Primary state-dependent model", data = "adj", incident = "I13",
                 protection = "R", timing = "adjacent", extras = character(0), weighted = FALSE),
  S1 = list(label = "S1 broader incident state (incident12)", data = "adj", incident = "I12",
            protection = "R", timing = "adjacent", extras = character(0), weighted = FALSE),
  S2 = list(label = "S2 monitoring-excluded protection (R3)", data = "adj", incident = "I13",
            protection = "R3", timing = "adjacent", extras = character(0), weighted = FALSE),
  S3 = list(label = "S3 next-wave web mode", data = "adj", incident = "I13",
            protection = "R", timing = "adjacent", extras = "web_mode_next", weighted = FALSE),
  S4 = list(label = "S4 strict Wave 3 to Wave 4 timing", data = "s4", incident = "I13",
            protection = "R", timing = "none", extras = character(0), weighted = FALSE),
  S5 = list(label = "S5 two-wave-gap timing", data = "gap2", incident = "I13",
            protection = "R", timing = "gap2", extras = character(0), weighted = FALSE),
  S6 = list(label = "S6 sector-adjusted (15 categories)", data = "adj", incident = "I13",
            protection = "R", timing = "adjacent", extras = "sector", weighted = FALSE),
  S7 = list(label = "S7 observable-retention weighted", data = "adj", incident = "I13",
            protection = "R", timing = "adjacent", extras = character(0), weighted = TRUE),
  S8 = list(label = "S8 broad-region-adjusted (6 categories)", data = "adj", incident = "I13",
            protection = "R", timing = "adjacent", extras = "region2", weighted = FALSE)
)
ST_STATES <- list(R = c(0.75, 1.00), R3 = c(2 / 3, 1))

st_parse_args <- function(argv) {
  get <- function(flag) {
    i <- match(flag, argv)
    if (is.na(i) || i == length(argv)) stop(paste("missing argument", flag))
    argv[i + 1]
  }
  list(mode = "reproduction", input_dir = normalizePath(get("--input-dir"), mustWork = TRUE),
       output_dir = normalizePath(get("--output-dir"), mustWork = TRUE))
}

st_sha256 <- function(path) digest::digest(file = path, algo = "sha256")

# The guard hashes the input files before any data file is parsed and writes nothing when a hash
# differs from the expected value.
st_guard <- function(args) {
  hashes <- vapply(names(ST_EXPECTED_INPUTS), function(f) st_sha256(file.path(args$input_dir, f)), "")
  if (!all(hashes == ST_EXPECTED_INPUTS)) stop("REFUSED: input SHA256 differs from the expected value")
  hashes
}

st_count_checks <- function(adj, gap2, elig, mode) {
  s4 <- adj$wave_t == 3
  chk <- list(
    adjacent_pairs = c(nrow(adj), 447), adjacent_organizations = c(length(unique(adj$id)), 308),
    adjacent_incident_free = c(sum(adj$I13_t == 0), 190),
    adjacent_incident_positive = c(sum(adj$I13_t == 1), 257),
    s4_pairs = c(sum(s4), 89), s4_incident_free = c(sum(s4 & adj$I13_t == 0), 35),
    gap2_pairs = c(nrow(gap2), 154), gap2_organizations = c(length(unique(gap2$id)), 116),
    gap2_incident_free = c(sum(gap2$I13_t == 0), 64),
    retention_rows = c(nrow(elig), 686), retention_retained = c(sum(elig$retain_adjacent), 447)
  )
  if (mode == "reproduction") {
    chk <- c(chk, list(
      adjacent_events = c(sum(adj$I13_next), 268),
      adjacent_onset_events = c(sum(adj$I13_next[adj$I13_t == 0]), 73),
      s4_onset_events = c(sum(adj$I13_next[s4 & adj$I13_t == 0]), 15),
      gap2_events = c(sum(gap2$I13_next), 105),
      gap2_onset_events = c(sum(gap2$I13_next[gap2$I13_t == 0]), 32)
    ))
  }
  out <- data.frame(check = names(chk), observed = vapply(chk, `[`, 0, 1),
                    expected = vapply(chk, `[`, 0, 2))
  out$pass <- out$observed == out$expected
  out
}

st_value_checks <- function(df, name) {
  bin <- intersect(c("I13_t", "I13_next", "I12_t", "I12_next", "web_mode_next"), names(df))
  ok <- !anyNA(df) && all(unlist(df[bin]) %in% c(0, 1)) &&
    all(df$R >= 0 & df$R <= 1) && all(df$R3 >= 0 & df$R3 <= 1) && all(df$E >= 0 & df$E <= 1)
  data.frame(check = paste0(name, "_values_valid"), observed = as.numeric(ok), expected = 1, pass = ok)
}

# Outcome-blind retention model used for the S7 weights.
st_retention <- function(elig, adj) {
  fr <- st_fit_glm(retain_adjacent ~ factor(wave) + I13 + E + R, elig)
  if (inherits(fr$fit, "st_fit_error") || !isTRUE(fr$fit$converged)) {
    return(list(ok = FALSE, reason = "retention model failed or did not converge"))
  }
  p <- unname(predict(fr$fit, newdata = data.frame(wave = adj$wave_t, I13 = adj$I13_t,
                                                   E = adj$E, R = adj$R), type = "response"))
  if (any(!is.finite(p)) || any(p <= 0) || any(p >= 1)) {
    return(list(ok = FALSE, reason = "invalid retention probability"))
  }
  w <- 1 / p
  list(ok = TRUE, fit = fr$fit, p = p, w = w / mean(w), warnings = fr$warnings)
}

st_run_spec <- function(key, spec, datasets, retention) {
  smp <- st_build_sample(datasets[[spec$data]], spec$incident, spec$protection, spec$timing, spec$extras)
  d <- smp$data
  base <- list(key = key, label = spec$label, sample = smp, n = nrow(d), G = length(unique(d$id)),
               n_it0 = sum(d$It == 0), n_it1 = sum(d$It == 1), events = sum(d$y),
               events_it0 = sum(d$y[d$It == 0]), states = ST_STATES[[spec$protection]],
               ok = FALSE, reasons = character(0), warnings = character(0), separation = FALSE)
  if (spec$weighted && !retention$ok) {
    base$reasons <- paste("retention weights unavailable:", retention$reason)
    return(base)
  }
  w <- if (spec$weighted) retention$w else NULL
  fr <- st_fit_glm(smp$formula, d, w)
  chk <- st_check_fit(fr, d$id)
  base$warnings <- fr$warnings
  base$separation <- fr$separation
  base$fit <- fr$fit
  if (!chk$ok) {
    base$reasons <- chk$reasons
    return(base)
  }
  avg_w <- if (spec$weighted) retention$w else rep(1, nrow(d))
  ct <- tryCatch(st_all_contrasts(fr$fit, chk$vc, d$It, base$states[1], base$states[2], smp$r_center, avg_w),
                 error = function(e) conditionMessage(e))
  if (is.character(ct)) {
    base$reasons <- paste("contrast failure:", ct)
    return(base)
  }
  base$ok <- TRUE
  base$vc <- chk$vc
  base$coefs <- st_coef_table(fr$fit, chk$vc)
  base$contrasts <- ct
  base
}

st_get <- function(res, term, field) res$coefs[res$coefs$term == term, field]
st_ct <- function(res, name, field) res$contrasts[res$contrasts$contrast == name, field]

st_decide <- function(results) {
  prim <- results$PRIMARY
  sens <- results[setdiff(names(results), "PRIMARY")]
  contra <- vapply(sens, function(r) r$ok && st_is_contradictory(st_get(r, "Rc", "lcl"), st_ct(r, "O", "estimate")), TRUE)
  dec <- list(contradictory_specs = names(contra)[contra],
              non_estimable_specs = names(sens)[!vapply(sens, `[[`, TRUE, "ok")],
              integrity_ok = isTRUE(prim$ok), conditions = NULL, influence_triggered = FALSE,
              influence_class = "NOT_TRIGGERED", influence_breaches = NA_integer_)
  if (!dec$integrity_ok) {
    dec$verdict <- ST_VERDICT_FAIL
    dec$failure_reason <- paste("primary specification failed:", paste(prim$reasons, collapse = "; "))
  } else {
    dec$conditions <- st_support_conditions(
      st_get(prim, "Rc", "estimate"), st_get(prim, "Rc", "ucl"),
      st_ct(prim, "O", "estimate"), st_ct(prim, "O", "lcl"),
      st_get(prim, "It:Rc", "estimate"), st_get(prim, "It:Rc", "lcl"), st_get(prim, "It:Rc", "ucl"),
      st_ct(prim, "S", "estimate"), st_ct(prim, "S", "lcl"), st_ct(prim, "S", "ucl"),
      any(contra))
    dec$verdict <- st_classify(dec$conditions, st_get(prim, "Rc", "estimate"), st_ct(prim, "O", "estimate"), TRUE)
    dec$failure_reason <- ""
  }
  dir_ok <- function(r) r$ok && st_get(r, "Rc", "estimate") < 0 && st_ct(r, "O", "estimate") > 0
  timing_clean <- !any(c("S4", "S5") %in% dec$contradictory_specs)
  dec$timing_onset_language_allowed <- dec$verdict != ST_VERDICT_FAIL && timing_clean &&
    (dir_ok(results$S4) || dir_ok(results$S5))
  same_sign <- function(r) {
    if (!r$ok || !prim$ok) return(FALSE)
    a <- st_get(prim, "It:Rc", "estimate")
    b <- st_get(r, "It:Rc", "estimate")
    abs(a) > ST_SIGN_TOL && abs(b) > ST_SIGN_TOL && sign(a) == sign(b)
  }
  dec$state_dependence_temporally_replicated <- dec$verdict == ST_VERDICT_SUPPORT && timing_clean &&
    (same_sign(results$S4) || same_sign(results$S5))
  dec
}

# An influence FAIL blocks a claim under PARTIAL as well as SUPPORT.
st_claim_status <- function(dec) {
  if (dec$verdict == ST_VERDICT_FAIL) return("NO_CLAIM")
  if (dec$verdict == ST_VERDICT_PARTIAL) {
    if (identical(dec$influence_class, "FAIL")) return("NO_CLAIM_INFLUENCE_FAIL")
    return("PARTIAL_DIRECTIONAL_ONSET_ONLY_NO_STATE_DEPENDENCE_CLAIM")
  }
  switch(dec$influence_class,
         PASS = "BOUNDED_SUPPORT_CLAIM",
         PARTIAL_TRIGGER = "SUPPORT_WITH_INFLUENCE_DISCLOSURE",
         "NO_CLAIM_INFLUENCE_FAIL")
}

st_main <- function() {
  args <- st_parse_args(commandArgs(trailingOnly = TRUE))
  hashes <- st_guard(args)
  started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  out <- function(name) file.path(args$output_dir, name)
  state <- list(args = args, hashes = hashes, started = started, results = list(), checks = NULL,
                retention = list(ok = FALSE, reason = "not computed"), decision = NULL, influence = NULL)
  status <- tryCatch({
    adj <- read.csv(file.path(args$input_dir, "configuration_adjacent_pairs.csv"), check.names = FALSE)
    gap2 <- read.csv(file.path(args$input_dir, "configuration_two_wave_gap_pairs.csv"), check.names = FALSE)
    elig <- read.csv(file.path(args$input_dir, "configuration_retention_eligibility.csv"), check.names = FALSE)
    state$checks <- rbind(st_count_checks(adj, gap2, elig, args$mode),
                          st_value_checks(adj, "adjacent"), st_value_checks(gap2, "gap2"))
    if (!all(state$checks$pass)) stop("pre-estimation count or value check failed")
    state$retention <- st_retention(elig, adj)
    datasets <- list(adj = adj, s4 = adj[adj$wave_t == 3, , drop = FALSE], gap2 = gap2)
    for (key in names(ST_SPECS)) {
      state$results[[key]] <- tryCatch(
        st_run_spec(key, ST_SPECS[[key]], datasets, state$retention),
        error = function(e) list(key = key, label = ST_SPECS[[key]]$label, ok = FALSE,
                                 reasons = paste("unexpected error:", conditionMessage(e)),
                                 warnings = character(0), separation = FALSE))
    }
    state$decision <- st_decide(state$results)
    if (state$decision$verdict %in% c(ST_VERDICT_SUPPORT, ST_VERDICT_PARTIAL)) {
      prim <- state$results$PRIMARY
      state$influence <- st_influence(prim$sample, st_get(prim, "Rc", "estimate"), st_get(prim, "It:Rc", "estimate"))
      nb <- sum(state$influence$breach)
      state$decision$influence_triggered <- TRUE
      state$decision$influence_breaches <- nb
      state$decision$influence_class <- st_influence_class(nb)
    }
    state$decision$claim_status <- st_claim_status(state$decision)
    "COMPLETED"
  }, error = function(e) {
    state$decision <<- list(verdict = ST_VERDICT_FAIL, integrity_ok = FALSE,
                            failure_reason = paste("execution integrity failure:", conditionMessage(e)),
                            claim_status = "NO_CLAIM", influence_class = "NOT_TRIGGERED",
                            influence_triggered = FALSE, influence_breaches = NA_integer_,
                            contradictory_specs = character(0), non_estimable_specs = character(0),
                            timing_onset_language_allowed = FALSE,
                            state_dependence_temporally_replicated = FALSE, conditions = NULL)
    "INTEGRITY_FAILURE"
  })
  state$status <- status
  state$finished <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  state$script_hashes <- vapply(c("state_transition_functions.R", "state_transition_outputs.R", "state_transition_models.R"),
                                function(f) st_sha256(file.path(st_script_dir, f)), "")
  wrote <- tryCatch({
    st_write_outputs(state, out)
    TRUE
  }, error = function(e) {
    st_write_fallback(state, out, conditionMessage(e))
    FALSE
  })
  message("State-transition run status: ", status, "; verdict: ", state$decision$verdict, "; outputs written: ", wrote)
  if (status != "COMPLETED" || !wrote) quit(status = 3)
}

# Run only as a script; sourcing defines the functions without executing.
if (sys.nframe() == 0L) st_main()

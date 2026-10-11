suppressPackageStartupMessages(library(survey))

input_csv <- "results/uk/build/uk_analysis_input.csv"
output_dir <- "results/uk/build"

d <- read.csv(input_csv, stringsAsFactors=FALSE)

required <- c(
  "uniqser","IMID2","year","year_calendar","weight","w_norm",
  "sizeb","size_medium","region_comb","sector_comb2",
  "online4","online5","online9","mobile","mobile_bin",
  "footprint_score","footprint_breadth",
  "rules2","rules4","rules5","rules7",
  "response_score","response_breadth",
  "manage1","manage2"
)

checks <- c(
  row_count = nrow(d) == 2859,
  required_columns = all(required %in% names(d)),
  primary_m0_complete = all(complete.cases(d[,c("response_breadth","footprint_breadth","year_calendar","size_medium","region_comb","w_norm")])),
  manage1_valid = all(d$manage1 %in% c(0,1)) && !any(is.na(d$manage1)),
  manage2_valid = all(d$manage2 %in% c(0,1)) && !any(is.na(d$manage2)),
  weights_positive = all(is.finite(d$w_norm)) && all(d$w_norm > 0),
  year_levels_valid = identical(sort(unique(d$year_calendar)), 2017:2021),
  size_levels_valid = identical(sort(unique(d$size_medium)), c(0L,1L)),
  region_complete = !any(is.na(d$region_comb)),
  footprint_breadth_range = all(d$footprint_breadth >= 0 & d$footprint_breadth <= 1),
  response_breadth_range = all(d$response_breadth >= 0 & d$response_breadth <= 1),
  footprint_score_range = all(d$footprint_score %in% 0:4),
  response_score_range = all(d$response_score %in% 0:4)
)

if (!all(checks)) {
  writeLines(paste(names(checks), checks, sep="="), file.path(output_dir, "uk_design_checks.txt"))
  quit(status=2)
}

d$year_factor <- factor(d$year_calendar, levels=2017:2021)
region_levels <- sort(unique(d$region_comb))
d$region_factor <- factor(d$region_comb, levels=region_levels)

design <- svydesign(ids=~1, weights=~w_norm, data=d)

checks <- c(
  checks,
  design_created = inherits(design, "survey.design2") || inherits(design, "survey.design"),
  design_rows = nrow(design$variables) == 2859,
  design_degf = degf(design) == 2858,
  year_reference_2017 = levels(d$year_factor)[1] == "2017",
  region_reference_smallest = as.numeric(levels(d$region_factor)[1]) == min(region_levels)
)

writeLines(paste(names(checks), checks, sep="="), file.path(output_dir, "uk_design_checks.txt"))
writeLines(c(
  paste("analysis_rows", nrow(d), sep="="),
  paste("design_rows", nrow(design$variables), sep="="),
  paste("design_degf", degf(design), sep="="),
  paste("year_levels", paste(levels(d$year_factor), collapse=","), sep="="),
  paste("year_reference", levels(d$year_factor)[1], sep="="),
  paste("region_levels", paste(levels(d$region_factor), collapse=","), sep="="),
  paste("region_reference", levels(d$region_factor)[1], sep="="),
  paste("size_values", paste(sort(unique(d$size_medium)), collapse=","), sep="="),
  paste("manage1_values", paste(sort(unique(d$manage1)), collapse=","), sep="="),
  paste("manage2_values", paste(sort(unique(d$manage2)), collapse=","), sep="="),
  paste("survey_version", as.character(packageVersion("survey")), sep="="),
  paste("r_version", as.character(getRversion()), sep="=")
), file.path(output_dir, "uk_design_manifest.txt"))

saveRDS(d, file.path(output_dir, "uk_analysis_ready.rds"), version=3)
saveRDS(design, file.path(output_dir, "uk_survey_design.rds"), version=3)

if (!all(checks)) {
  quit(status=3)
}

cat("UK_DESIGN_CHECKS_PASS\n")
cat("analysis_rows=", nrow(d), "\n", sep="")
cat("design_degf=", degf(design), "\n", sep="")
cat("year_reference=", levels(d$year_factor)[1], "\n", sep="")
cat("region_reference=", levels(d$region_factor)[1], "\n", sep="")

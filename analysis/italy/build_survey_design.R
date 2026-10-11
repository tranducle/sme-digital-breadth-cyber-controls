suppressPackageStartupMessages(library(survey))

expected_r <- "4.6.1"
expected_survey <- "4.5"

if (as.character(getRversion()) != expected_r) warning("Results were produced with R 4.6.1")
if (as.character(packageVersion("survey")) != expected_survey) warning("Results were produced with survey 4.5")

d <- read.csv("results/italy/build/italy_analysis_input.csv", stringsAsFactors=FALSE)

required <- c(
  "codice_","size_medium","rip","mac","dom4","weight",
  "footprint_score","footprint_breadth",
  "response_score","response_breadth",
  "digital_system_breadth_6","security_core_breadth_6"
)

checks <- c(
  row_count = nrow(d) == 16809,
  required_columns = all(required %in% names(d)),
  primary_complete = all(complete.cases(d[,c("response_breadth","footprint_breadth","size_medium","rip","mac","weight")])),
  positive_weights = all(is.finite(d$weight)) && all(d$weight > 0),
  footprint_range = all(d$footprint_breadth >= 0 & d$footprint_breadth <= 1),
  response_range = all(d$response_breadth >= 0 & d$response_breadth <= 1),
  footprint_score_range = all(d$footprint_score %in% 0:4),
  response_score_range = all(d$response_score %in% 0:10),
  size_values = identical(sort(unique(d$size_medium)), c(0L,1L))
)

if (!all(checks)) {
  writeLines(paste(names(checks),checks,sep="="), "results/italy/build/italy_design_checks.txt")
  quit(status=2)
}

d$size_factor <- factor(d$size_medium, levels=c(0,1), labels=c("Small","Medium"))
d$region_factor <- factor(d$rip, levels=sort(unique(d$rip)))
d$macrosector_factor <- factor(d$mac, levels=sort(unique(d$mac)))
d$ict_domain_factor <- factor(d$dom4, levels=sort(unique(d$dom4)))

design <- svydesign(ids=~1, weights=~weight, data=d)

checks <- c(
  checks,
  design_created = inherits(design,"survey.design2") || inherits(design,"survey.design"),
  design_rows = nrow(design$variables) == 16809,
  design_degf = degf(design) == 16808,
  size_reference = levels(d$size_factor)[1] == "Small",
  region_reference = levels(d$region_factor)[1] == "ITC",
  macrosector_reference = levels(d$macrosector_factor)[1] == "M1"
)

writeLines(paste(names(checks),checks,sep="="), "results/italy/build/italy_design_checks.txt")
writeLines(c(
  paste("analysis_rows",nrow(d),sep="="),
  paste("design_rows",nrow(design$variables),sep="="),
  paste("design_degf",degf(design),sep="="),
  paste("size_levels",paste(levels(d$size_factor),collapse=","),sep="="),
  paste("region_levels",paste(levels(d$region_factor),collapse=","),sep="="),
  paste("macrosector_levels",paste(levels(d$macrosector_factor),collapse=","),sep="="),
  paste("ict_domain_levels",paste(levels(d$ict_domain_factor),collapse=","),sep="="),
  paste("r_version",as.character(getRversion()),sep="="),
  paste("survey_version",as.character(packageVersion("survey")),sep="=")
), "results/italy/build/italy_design_manifest.txt")

saveRDS(d, "results/italy/build/italy_analysis_ready.rds", version=3)
saveRDS(design, "results/italy/build/italy_survey_design.rds", version=3)

if (!all(checks)) quit(status=3)

cat("ITALY_DESIGN_CHECKS_PASS\n")
cat("n=",nrow(d),"\n",sep="")
cat("degf=",degf(design),"\n",sep="")
cat("region_reference=",levels(d$region_factor)[1],"\n",sep="")
cat("macrosector_reference=",levels(d$macrosector_factor)[1],"\n",sep="")

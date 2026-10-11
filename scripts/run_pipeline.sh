#!/usr/bin/env bash
# Run the full analysis pipeline from the repository root.
# Raw data must already be in place under data/raw/ (see DATA.md).
set -euo pipefail
cd "$(dirname "$0")/.."

PY=${PYTHON:-python3}
RS=${RSCRIPT:-Rscript}

mkdir -p results/uk/{build,primary,sensitivity,decomposition} \
         results/italy/{build,models} results/size_gap results/measurement \
         results/panel results/next_wave/{inputs,configuration,state_transition} \
         results/table_inputs results/figure_inputs

step() { echo "==> $*"; }

step "United Kingdom: analysis dataset and survey design"
$PY analysis/uk/build_analysis_dataset.py
$RS analysis/uk/build_survey_design.R

step "United Kingdom: primary models"
$RS analysis/uk/model_m0_baseline.R
$RS analysis/uk/model_m1_primary.R
$RS analysis/uk/model_m2_management_interactions.R
$RS analysis/uk/model_m3_holm_correction.R

step "United Kingdom: sensitivity analyses"
$PY analysis/uk/sensitivity/build_unweighted_input.py
for s in s1_unweighted s2_original_weights s3_sector_adjusted s4_categorical_breadth \
         s5_linear_count s6_joint_organization s7_leave_one_wave_out d1_later_window_without_sector; do
  $RS "analysis/uk/sensitivity/${s}.R"
done
$PY analysis/uk/sensitivity/d2_attenuation_comparison.py

step "Italy: analysis dataset, survey design and primary model"
$PY analysis/italy/build_analysis_dataset.py
$RS analysis/italy/build_survey_design.R
$RS analysis/italy/model_primary.R

step "Italy: sensitivity analyses"
for s in s1_equal_weights s2_categorical_breadth s3_broader_technology s4_core_response \
         s5_linear_response s6_ict_sector_adjusted s7_high_weight_influence; do
  $RS "analysis/italy/sensitivity/${s}.R"
done
$PY analysis/italy/sensitivity/build_all_sizes_input.py
$RS analysis/italy/sensitivity/s8_all_firm_sizes.R

step "Firm-size gap diagnostics"
$RS analysis/size_gap/size_gap_specifications.R
$RS analysis/size_gap/size_gap_attenuation.R

step "Measurement boundary check"
$PY analysis/measurement/alignment_audit.py

step "Longitudinal panel"
$PY analysis/panel/export_panel_columns.py
$RS analysis/panel/panel_models.R

step "Next-wave analyses"
$PY analysis/next_wave/build_inputs.py
$RS analysis/next_wave/configuration_models.R
$RS analysis/next_wave/state_transition_models.R \
    --input-dir results/next_wave/inputs --output-dir results/next_wave/state_transition

step "Inputs for manuscript tables and figures"
$RS analysis/reporting/table_inputs.R
$RS analysis/reporting/figure_inputs.R

echo "Pipeline finished."

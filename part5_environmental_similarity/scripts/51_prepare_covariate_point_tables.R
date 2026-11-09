# Part 2: reconcile and prepare the downloaded GEE point-extraction tables.
# This step does not fit a model and does not modify the raw CSV exports.

source(file.path("R", "covariates", "prepare_point_tables.R"))

covariate_point_tables <- prepare_covariate_point_tables(getwd())
print(covariate_point_tables$availability_summary)
message(
  "Prepared covariate point tables under data/derived/covariate_analysis/. ",
  "Review outputs/covariate_analysis/part2_qa/POINT_TABLE_SUMMARY.md."
)

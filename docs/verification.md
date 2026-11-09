# Organization checks

Checked on 30 October 2026 using R 4.6.0 and the supplied six-core simulated workbook.

- All five numbered scripts pass R syntax parsing.
- `script_1.R` completed from the project folder. It produced one mapping table with six core records, with 6/5/4 usable stocks in the 0–30/0–50/0–100 reporting groups. The profile charts, comparison figures and Part 3 summary were generated under `outputs/script1/`.
- `script_2.R` completed in a separate R process using Script 1's generated mapping table. Area polygons, the boundary, four IDW rasters, tuning tables, totals and the combined map were written under `outputs/script2/`.
- Script 3's file-reading section was evaluated in a fresh environment. It loaded six cores and the study boundary from the canonical file paths. Model fitting was not run.
- Script 4's file-reading section was evaluated in a fresh environment. It loaded the current field table, area polygons, study boundary and IDW maps from the new paths. Its full comparison calculations were not run.
- Script 5 was syntax-checked; its full analysis was not run.
- Reference preparation was rebuilt in a separate copy using the renamed source files and portable tools. The decompressed reference values, eligibility records, ecosystem CSV and all five summary tables matched the supplied prepared results.
- Copied source data retained their recorded checksums. All 87 source files, excluding macOS `.DS_Store` files, were preserved unchanged under `archive/transfer_original/`.
- The active scripts contain no `Draft1` or Windows `Y:/` paths. No CSVs or rasters are loose in the project folder.

These checks validate the organization and the stated execution paths. They do not establish scientific validity for model or uncertainty methods that are still to be completed in later workshop work.

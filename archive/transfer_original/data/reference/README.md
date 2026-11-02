# Script 1 carbon reference data

Prepared locally on 30 October 2026 from the five soil-core collections already used in the Carbon State of Data comparisons: NRCan peat, WoSIS Canada, NPDB, Janousek coastal cores, and Riley Ontario peat records. No additional catalogued datasets were added. This folder has not been uploaded to GitHub.

## Depth groups agreed with the user

| Group | Profile entry rule | Material included |
| --- | --- | --- |
| 0–30 cm | Sampled bottom reaches or exceeds 30 cm | Measured material from 0 to 30 cm |
| 0–50 cm | Sampled bottom is deeper than 30 cm | Measured material from 0 to the shallower of 50 cm and sampled bottom |
| 0–100 cm | Sampled bottom is deeper than 50 cm | Measured material from 0 to the shallower of 100 cm and sampled bottom |
| total | Any sampled depth | Material between recorded sampled top and bottom, where the required measurements are usable |

A 45 cm profile contributes its observed 0–45 cm stock to the 0–50 group. A 60 cm profile contributes its observed 0–60 cm stock to the 0–100 group. These are capped-stock groups and include partial target-depth coverage. They are not all complete observations to 50 or 100 cm. Each record retains its actual sampled depth and `complete_target_depth` flag. Summary tables record complete and partial target counts separately. Total sampled stocks span different depths and do not establish whole-soil-profile stocks or mineral contact.

For each section:

`stock contribution = section stock / section thickness × centimetres inside reporting window`

For example, a 30–60 cm section with 12 kg C/m² contributes 12/30 × 20 = 8 kg C/m² to the 0–50 group, added to the measured stock above 30 cm. This assumes a uniform carbon stock per centimetre within each measured section. There is no extrapolation below the sampled bottom. Midpoints can locate observations on profile plots; section boundaries determine the contribution to a stock sum.

## Files to read in R

- `reference_depth_values.csv.gz`: individual eligible profile values for the four depth groups. Variables are `stock_kgm2` (kg C/m²), `soc_gkg` (g C/kg) and `bd` (g/cm³). Concentration and density are thickness-weighted over the included measured interval.
- `summary_by_dataset.csv`, `summary_by_ecosystem.csv`, `summary_by_general_region.csv`, `summary_by_source_region.csv`, `summary_by_study_group.csv`: mean, median, SD, interquartile range, range, profile-record counts, coordinate-pair counts, actual depth ranges, and full/partial target counts.
- `depth_group_counts.csv`: eligible counts for every collection, group and variable, including zero counts.
- `profile_eligibility.csv.gz`: all profile/group/variable eligibility decisions and reasons. Every source profile is represented, including those without computable carbon stock.
- `resources.csv`: exactly five resource entries, source links, availability notes, original terms and preparation cautions.
- `profile_metadata.csv.gz`: existing profile metadata, original collection/study groupings, ecosystem and region labels.
- `cores/`: the prepared source layer/profile files and supporting NRCan/Janousek metadata. Compression preserves each original CSV byte for byte. R's `readr::read_csv()` reads `.csv.gz` directly.

Example:

```r
reference <- readr::read_csv("data/reference/reference_depth_values.csv.gz")
resources <- readr::read_csv("data/reference/resources.csv")
```

## Preparation and interpretation

The build uses the same source-stock fields and conversions as the previous SOD interval method, recorded under `provenance/previous_interval_method.R`. NRCan provider sample stocks and NPDB stocks are converted from Mg C/ha to kg C/m² by dividing by 10. WoSIS and Janousek prepared stocks are already kg C/m². NRCan SOC percentages and NPDB/Janousek/Riley SOC percentages are converted to g C/kg by multiplying by 10. Prepared WoSIS organic carbon is already g C/kg.

Only the surface-datum comparison is generated here. The inherited convention rebases NPDB's organic-inclusive profiles to their recorded top; other collections retain their supplied depth coordinates. Surface and mineral-datum results are not pooled. Riley bulk-density units remain unverified, so it contributes organic-carbon concentration only, with no computed density or stock values. Source stocks are not reconstructed by filling missing density. NRCan measured/estimated flags remain in the source layers. Janousek profiles with inconsistent source thickness/boundaries are excluded as in the previous preparation.

Invalid or zero-width source layers are excluded from integration and counted in the preparation log. Gaps and overlaps in a reporting window prevent a quantitative value for that window. This build uses a floating-point tolerance of 1e-8 cm rather than the old 0.5 cm coverage allowance. Negative stock/concentration and nonpositive density are excluded. These rules and the new user-selected eligibility cutoffs mean the counts need not reproduce the earlier summary snapshot. All original prepared source CSVs remain available separately.

Counts refer to profile records, not necessarily independent sites. Cross-collection duplication has not been resolved. Missing ecosystem classifications and unlinked original study attribution are retained explicitly. Janousek spans Pacific coastal settings from Alaska to Mexico; the five collections are broad descriptive context rather than interchangeable local analogues or independent validation of the simulated community data.

To rebuild from these local files:

```text
python3 data/reference/provenance/build_depth_groups.py
```

`provenance/source_manifest.csv` records source and stored checksums. `provenance/depth_group_settings.json` records the depth rules. `provenance/build_log.txt` records counts. The paste-ready Script 1 stock calculation and comparison were checked on the six-core example, including a section crossing 50 cm, exact cutoffs, incomplete cores, and missing/gapped/overlapping sections.

## Source terms and public repository inclusion

This is a local staging bundle and does not assign a common licence to third-party data. NRCan and NPDB provider records identify Open Government Licence – Canada. Prepared WoSIS rows retain their original licence fields, including public domain, CC BY and CC BY-NC records. [ISRIC's WoSIS FAQ](https://docs.isric.org/globaldata/wosis/faq-wosis.html) confirms that contributing-provider conditions apply. The [Janousek archive](https://doi.org/10.25573/serc.28127486.v1) needs its original attribution and reuse terms retained. Riley's exact original release and redistribution terms remain unresolved. The full local folder should not be described as a cleared unrestricted public release.

For the intended GitHub layout, this folder belongs at `Transfer_OCT_29/data/reference/` on the `Transfer_Oct` branch. Script 1 is then run with `Transfer_OCT_29` as the working directory. Public upload was not performed during this task.

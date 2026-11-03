# Reference observations

This folder contains the five soil-core collections used in the Carbon State of Data comparisons: NRCan peat, WoSIS Canada, NPDB, Janousek coastal cores and Riley Ontario records.

## Main files

- `reference_depth_values.csv.gz`: individual usable stock, carbon concentration and bulk-density values for the 0–30, 0–50, 0–100 and total reporting groups.
- `ecosystem_carbon_stocks.csv`: the single stock-only CSV, including ecosystem, source collection, source URL, study group, actual sampled depths and partial-depth flags.
- `summary_by_*.csv`: descriptive summaries by collection, ecosystem, general region, source region and study group.
- `resources.csv`: source links, availability, terms and preparation notes for the five collections.
- `depth_group_counts.csv`: usable record counts, including zero counts.
- `profile_eligibility.csv.gz`: each profile's inclusion/exclusion decisions and reasons.
- `profile_metadata.csv.gz`: the recorded ecosystem, region and study metadata.
- `cores/`: original prepared layer/profile CSVs, compressed without changing their values. Provider-specific IDs and recorded licence fields are retained.
- `provenance/`: source checksums, preparation settings and historical source-method notes. Files beginning `previous_` describe the earlier preparation and can contain its original paths and file names.

R's `readr::read_csv()` opens `.csv.gz` directly. Users do not need to unzip these tables.

## Methods retained

The reporting groups retain the agreed entry depths: reach 30 for 0–30; exceed 30 for 0–50; exceed 50 for 0–100; retain sampled profile material for total. Stock is allocated by section overlap, assuming uniform carbon stock per centimetre within each section. Concentration and density are thickness-weighted. Missing measurements are not filled and unsampled depths are not extrapolated. Actual sampled depths and complete/partial target flags remain separate from the reporting-group label.

Source stock fields and unit conversions are retained from the previous SOD preparation. NRCan/NPDB stock in Mg C/ha is divided by 10 to give kg C/m². Prepared WoSIS and Janousek stocks are already kg C/m². Carbon percentages are multiplied by 10 for g C/kg. Riley density units remain unverified, so that collection contributes concentration but no computed density or stock.

The comparison uses the surface datum. NPDB's organic-inclusive profiles retain the previous rebasing to their recorded profile top. Other collections retain their supplied depths. Surface and mineral results are not pooled. Invalid section depths are recorded, and gaps/overlaps prevent a quantitative value for the affected interval. Source boundary/thickness inconsistencies and measurement flags remain documented in provenance.

Counts are profile records and are not asserted to be independent sites. Cross-collection duplication is unresolved. Missing ecosystem labels remain unknown. Janousek covers coastal environments from Alaska to Mexico and is broader context for an inland site. Community percentile summaries pool the available reference records and are not an ecosystem-matched habitat assessment.

## Rebuild

From the project folder, run `python3 tools/prepare_reference.py`, then `python3 tools/export_ecosystem_stocks.py`. This uses the bundled source tables without a new download.

## Source terms

These are local prepared snapshots, not a new blanket data licence. NRCan and NPDB provider metadata identify Open Government Licence – Canada. WoSIS retains mixed provider-level public-domain, CC BY and CC BY-NC terms. Janousek retains the archive route and contributing-study information. Riley's exact original release and redistribution terms remain unresolved. Preserve the records in `resources.csv` and the source files when considering any later public distribution.

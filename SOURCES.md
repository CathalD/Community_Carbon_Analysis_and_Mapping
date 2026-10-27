# Data sources and comparison snapshot

## Community workbook

`run_example.R` imports the supplied workbook and recalculates its layer and core values. The simple `field_data_template.csv` contains the same 29 slices with the confirmed example density correction already applied. `run_my_data.R` reads a user's eight-column CSV (or first-sheet Excel table) and polygon boundary; the example correction is disabled for that route. Both preserve the input file and save derived CSVs separately. Spreadsheet macros are not required or executed.

The optional full course uses the workbook named in `config.R`. `course_config$data_origin` records how those observations should be described. Revised illustrative measurements and conceptual locations are examples; they do not establish actual ecosystem conditions at those coordinates.

Community stocks use dry bulk density in g/cm³, organic carbon entered as a percentage, and interval thickness in cm:

`stock_kg_m2 = bulk_density * organic_carbon_percent * thickness * 0.1`

The concentration/density comparison uses thickness-weighted interval averages. Carbon percentages are multiplied by 10 to compare with reference concentrations in g C/kg. Intervals must be completely sampled, with no overlaps or gaps; missing deeper measurements are not extrapolated. The default depth reference is the profile surface. Changing `depth_datum` requires changing the input depths to that datum, not just changing the label.

## Precomputed Carbon State of Data reference summaries

`data/reference/summary_by_*.csv.gz` contains five aggregate tables prepared by the existing Carbon SOD/local workshop workflow on 26 October 2026. The dimensions are collection, study or contributing collection, recorded ecosystem, Canadian ecozone/general region, and province/state/source region. The snapshot manifest records the original table names, row counts and checksums. No profile coordinates, individual reference measurements or raw Riley/WoSIS/NPDB metadata are included in this course bundle.

The aggregate tables cover five loaded collections. These are the eligible **surface-datum** profile-record counts for the three overview quantities; each other depth and datum has its own counts in the supplied tables.

| Collection | Stock 0–30 cm | Organic C concentration 0–5 cm | Bulk density 0–5 cm |
|---|---:|---:|---:|
| NRCan peat | 1,211 | 1,213 | 1,213 |
| WoSIS Canada | 14 | 7,975 | 21 |
| National Pedon Database (NPDB) | 338 | 5,650 | 539 |
| Janousek Pacific coastal collection | 1,230 | 1,240 | 1,240 |
| Riley Ontario preparation | — | 279 | — |

Counts are eligible profile records, not necessarily independent sites. The tables also retain recorded coordinate-pair counts as an aggregate diagnostic; repeated locations and cross-collection duplicates have not been resolved. Small groups remain in the complete tables. Other catalogued data resources without loaded observations are outside this snapshot.

### What this download recalculates

`scripts/02_compare.R` reads the aggregate snapshot, recalculates the community quantities from its imported layers, and matches quantity, units, depth interval and depth datum. It then regenerates the comparison tables and figures. It **does not rebuild the reference statistics from raw reference profiles**. Raw-profile quality checks, alternative depth integration, reclassification, deduplication, resampling and hypothesis tests cannot be performed from these aggregates.

The original local workshop's `07_compare_ground_profiles.R` and its helpers reproduce the full reference summaries where the separately held raw SOD preparations are available. Those local inputs are not required for the downloadable MVP and are not copied into this repository. To refresh the public snapshot, rerun that local preparation, review the five summary tables, replace the compressed aggregate files and update the manifest.

### Reference preparation retained in the snapshot

- Concentration and density are overlap-weighted averages of the prepared source layers. The supported windows are 0–5, 5–15, 15–30, 30–60 and 60–100 cm.
- Stocks use prepared source layer-stock fields allocated by overlap for 0–30 and 0–100 cm. This preserves the SOD stock preparation rather than recalculating every reference from rounded concentration and density. Source stock and coarse-fragment conventions vary.
- The inherited reference coverage check allows a 0.5 cm tolerance and excludes gaps/overlaps outside that tolerance. It does not extrapolate beyond sampled coverage or fill missing density. Community interval checks are stricter.
- Surface and mineral datums remain separate. For NPDB, the surface preparation shifts the top of the profile to zero, including organic surface layers. The mineral preparation removes the organic layers and shifts the first remaining mineral horizon to zero. WoSIS and NPDB have mineral-datum reference summaries; the course's default community inputs use surface depths.
- Riley contributes organic carbon concentration only: its prepared bulk-density units remain unverified, so density and stock are excluded.
- Janousek profiles with recorded sample thickness inconsistent with the section boundaries were excluded. The collection covers Pacific coastal settings from Alaska to Mexico; it is broader context for an inland Ontario example.
- NRCan peat includes both measured and estimated properties. Numeric completeness does not mean that every underlying property was measured.

### Group labels and figure meaning

NRCan has recorded wetland classes, and Janousek has recorded coastal ecosystem codes. WoSIS and NPDB have no reliable ecosystem labels in the current preparation; `Not recorded` is preserved. Riley is labelled peatland without a verified subtype. Soil class is not converted into a forest/ecosystem label.

Janousek study groups are archive study codes, with paper citations not yet linked. WoSIS groups are contributing-collection codes, not confirmed individual papers. Original NRCan, NPDB and Riley study attribution is unlinked in this snapshot. Group labels state these limits.

Figure segments show the reference interquartile range; filled dots show medians and open diamonds show means. Dashed coloured lines show individual community cores at the same depth and datum. There are no reconstructed boxplot whiskers or raw reference points. The study and source-region stock figures show the 20 largest groups by record count; complete output tables retain every group. Comparisons against other ecosystem/region groups provide context, not a claim that those groups are local analogues or statistically different populations. They are not independent validation of the illustrative inputs.

## Original source routes and terms

| Collection | Source route | Important source facts |
|---|---|---|
| NRCan Peat Profile Database | [Open Canada provider record](https://open.canada.ca/data/en/dataset/e287b0dd-8b07-435b-8a76-0c7ec84ac228) | Provider record lists Open Government Licence – Canada. The historical local name `CanPeat` refers to this database, not the Waterloo Can-Peat portal. Original source references and measured/estimated flags matter. |
| WoSIS | [ISRIC WoSIS documentation](https://docs.isric.org/globaldata/wosis/) | Contributing-source terms vary. The snapshot does not grant a new licence to raw contributing observations. |
| National Pedon Database | [CanSIS NPDB](https://sis.agr.gc.ca/cansis/nsdb/npdb/index.html) | Prepared horizon and layer-stock conventions require source-specific interpretation. No raw NPDB release is bundled here. |
| Janousek coastal compilation | [Smithsonian archive, version 1](https://doi.org/10.25573/serc.28127486.v1) | Preserve archive attribution and contributing-study codes; not all records are Canadian. Consult the archive's original terms for raw data reuse. |
| Riley Ontario preparation | Exact original OGS release not resolved in the current preparation | Original source identification, density units and redistribution terms remain unresolved. No raw Riley observations are bundled. |

The course's code licence does not replace third-party data terms. Public accessibility of an archive, catalogue or hosted asset does not itself assign a common licence to all underlying sources. These aggregate comparison products retain their collection labels and the documented limitations above.

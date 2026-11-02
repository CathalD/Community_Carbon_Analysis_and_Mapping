# Carbon SOD ground-data snapshots

These local snapshots extend the workshop's optional ground-data comparison. They were copied from the existing Carbon_Data_SOD project and CommunityCarbonMap_v3 preparations on 26 October 2026. No values were changed during copying. Running Step 07 uses these files locally; it does not download data or require Earth Engine access.

## Files and roles

- `profile_intervals.csv`: SOD's 151,688 eligible depth-standardized values. Variables are organic C concentration in g/kg, dry bulk density in g/cm3 and stock in kg C/m2. Values are retained separately by depth interval and datum.
- `profile_locations.csv`: 20,690 profile IDs, recorded coordinates and Canadian ecozones.
- `profile_dataset_summary.csv`: the existing SOD preparation/coverage audit.
- `source_interval_method.R`: the original SOD interval-building method, retained as provenance. It references the original source-project paths and is not run by this workshop.
- `CANPEAT_CORES_DATA.csv`: original NRCan core metadata, including wetland class, province, sampling year, location notes and whole-profile measured percentages.
- `Global_Core_Locations.csv`: original Janousek archive core IDs, ecosystem codes, state/region codes and contributing study codes.
- `wosis_canada_layers_v2.csv`: Canadian WoSIS preparation with contributing dataset IDs and original method/licence fields.
- `npdb_canada_layers_v1_upload_fixed.csv`: NPDB preparation with province, year and soil-class fields. Soil class and drainage are not interpreted as ecosystem or forest-cover labels.
- `janousek_profiles_v1.csv`: prepared Janousek profile metadata, retained for source reconciliation.
- `Ground_data.csv`, `Studies.csv`, `Study_ground.csv`: the existing SOD catalogue and study-source links. A study's use of a database is not a profile-to-original-paper attribution and must not duplicate those observations.
- `NRCan_REFERENCES.csv`: provider publication references, downloaded from the public NRCan peat-database directory. The accompanying `SOURCEDATA.csv` core-to-reference crosswalk was not retrieved successfully; the code leaves original NRCan study attribution unlinked unless that file is later supplied.

Original SOD files reside in `/Users/cathaldoherty/Desktop/Carbon_Resources/Carbon_Data_SOD/Comparison/data/derived/` and `data/downloads/`. Layer preparations reside in the original project's `outputs/carbon_state_of_data/soil_uploads/`. Original NRCan and Janousek tables reside under `/Users/cathaldoherty/Desktop/Carbon_Resources/SoilDatabases_RStudio/SoilCarbonProfiles-main/`.

The comparison writes `input_manifest.csv` with file checksums. These are local working snapshots; original source-specific reuse terms continue to apply. Retain attribution and original licence fields when preparing a public workshop subset.

## Harmonization inherited from SOD

SOD uses overlap-weighted concentration and density, and sums each provider's layer stock field allocated by depth overlap. It does not fill missing density or extrapolate incomplete profiles. Its coverage checks use a 0.5 cm tolerance; the original method is retained above. This differs slightly from the stricter synthetic-core check and from deriving NRCan stocks afresh using rounded concentration/density fields in the local comparison.

- `surface`: depth from profile top, including organic layers. NPDB profiles are shifted from their original minimum depth to that profile top.
- `mineral`: organic surface horizons are removed from NPDB/WoSIS, then the first remaining mineral horizon is set to zero.
- Concentration/density windows: 0–5, 5–15, 15–30, 30–60 and 60–100 cm.
- Stock windows: 0–30 and 0–100 cm.
- Riley density and stock are excluded because density units remain unverified. Its organic carbon concentration remains available.
- Janousek profiles with inconsistent thickness versus section boundaries are excluded by the existing SOD method. The collection covers coastal sites across Alaska to Mexico and is not an exclusively Canadian or inland-forest dataset.
- Surface and mineral datums are never pooled. Synthetic overlays currently match only the surface datum.

## Metadata limits

NRCan provides wetland types; Janousek provides coastal ecosystem codes. WoSIS and NPDB preparations do not provide a reliable ecosystem classification, so their ecosystem is `Not recorded`. Riley is recorded only as peatland without a confirmed subtype.

Janousek's study IDs are actual source-archive study codes, with paper titles not yet joined. WoSIS dataset IDs are contributing-collection codes, not individual study publications. NRCan, NPDB and Riley original study attribution remains unlinked where no verified profile crosswalk is available. These distinctions are preserved in `study_attribution.csv`; no source names are inferred from coordinates or sampling-strategy codes.

General regions use the recorded Canadian ecozones. Outside/unmatched locations remain a separate group. Province/state spellings are normalized in `source_region` (for example ON and Ontario); original codes remain in `source_region_raw`. Janousek state codes are decoded within their country context, so CA means California rather than Canada. Distinct profile IDs and rounded coordinate pairs are not necessarily independent sampling sites. Cross-collection duplication is not resolved, so the outputs remain descriptive.

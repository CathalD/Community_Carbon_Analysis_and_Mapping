# Analysis-ready covariate point tables

Generated: 2026-09-24 13:17:05 UTC

- Raw Earth Engine exports were read without modification.
- Missing sentinel values (-9999) were converted to NA.
- External export: 1,564 unique Riley/CanPEAT cores.
- Complete external 0-30 cm target: 1,512 cores.
- Community export: 8 unique cores, all joined to Route L organic-depth contacts.
- All files contain the expected 119 predictors and five QA fields.

## Predictor availability

- external_all: 1564 rows; 1546 complete surface-context; 1563 complete AlphaEarth; 1556 complete terrain; 1538 complete across all 119 predictors.
- external_complete_0_30_target: 1512 rows; 1494 complete surface-context; 1511 complete AlphaEarth; 1504 complete terrain; 1486 complete across all 119 predictors.
- community: 8 rows; 5 complete surface-context; 8 complete AlphaEarth; 8 complete terrain; 5 complete across all 119 predictors.

Missing predictors remain NA and must be handled explicitly by each candidate model.

# Example workbook

`community_cores.xlsm` is the supplied revised illustrative workbook, preserved unchanged. R reads its input cells; it does not run macros or use embedded figures. No VBA project was present in this supplied file.

The author confirmed a correction to `3. Peat Data!J7`: **0.9 → 0.09 g/cm³**. `config.R` applies it on import and the QC output records it. Entered organic carbon is authoritative; conflicting LOI, water content and chemistry notes remain flagged and are not used to derive stock. The workbook's existing images are stale. Use the regenerated figures in `outputs/`.

All four core locations and ecosystem assignments remain hypothetical. These are teaching data, not a field survey or Fort Severn analysis.

For your own data, save a workbook in this format and change `workbook_path` in `config.R`. Set `density_override_J7 = NA_real_` to disable the example-specific correction. Update the core/plot IDs, coordinates, ecosystem labels and depth reference. Run Part 1 and review `qc_flags.csv` before mapping. The spatial modules deliberately require illustrative mode until you replace their simulated inputs and uncertainty assumptions.

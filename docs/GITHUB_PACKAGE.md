# Lightweight FS test GitHub snapshot

This branch contains the current four-method FS workflow, original workshop source archive, documentation, measured FS field tables/boundary, predictor-at-core tables, compact analysis/validation summaries, and current full-suite PNG figures. The analysis code and scientific assumptions are preserved. Notebook saved outputs were cleared without changing its code or markdown.

Large local inputs and regenerated artifacts are intentionally omitted: prior and covariate TIFFs, VRTs, GIS/model binaries, compressed reference profile tables, the large ecosystem reference table, per-cell prediction CSVs, execution/session caches, and duplicate PDFs. Their paths, sizes and source SHA-256 hashes are recorded in `github_package_manifest.csv`. The local originals were not deleted or changed.

## Reproduce the analysis

1. Open `FS_test.Rproj` and use an R library with the packages listed in the main README.
2. Restore omitted files to the same relative paths from the local FS test folder, using the manifest to check the original hashes. Full re-execution needs the omitted input data; reading the included summaries and figures does not.
3. If rebuilding cropped priors/covariates, edit `config.R`: `source_root` must point to the original CommunityCarbonMap_v3 directory containing `data/external`, `sothe_root` to the original Sothe release directory, and `covariate_root` to the completed downloaded covariate tiles. Existing defaults describe the source workstation; they do not automatically follow a moved folder. `prepare_inputs.R` rebuilds the local raster products. Reference preparation also needs the omitted source tables under `data/reference/cores`.
4. Run the model-stage commands in the main README after restoring/preparing inputs.

The included verification tables document the completed local run. They do not certify a data-incomplete checkout as runnable, or validate the ecological assumptions. `outputs/verification/input_manifest.csv` and `code_manifest.csv` describe the original analytical run; this GitHub package manifest separately records packaging omissions and the output-cleared notebook.

The four Li Gaussian-scale updates retain eight input records, including classified zeros, and remain sensitivity analyses. The current results include flagged negative spatial estimates and prior-retained regression cells. The README and run-results document retain these limitations. Source attribution and provider-specific terms are preserved in `data/reference/README.md` and `resources.csv`; this snapshot does not introduce a new blanket licence for third-party data.

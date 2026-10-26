# Reference summary inputs

The five compressed CSV tables are precomputed aggregate SOD reference statistics. R reads them directly using `gzfile()`; there is no separate decompression step or raw-data download.

`summary_snapshot_manifest.csv` records source filenames, snapshot date, row counts and checksums. All original variables, supported depth intervals and depth datums are retained.

See [SOURCES.md](../../SOURCES.md) for provenance, availability counts, source terms, interval methods and the distinction between regenerating the community comparisons and rebuilding the raw-reference summaries.
